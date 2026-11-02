// Interaction sol-radier par elements finis -- implementation de reference.
//
// Plaque   : element rectangulaire de Kirchhoff-Love ACM/MZC (4 noeuds, 12 ddl,
//            ddl nodaux (w, w_x, w_y)), coordonnees naturelles xi, eta dans [-1,1].
// Sol      : couches elastiques, souplesse obtenue en integrant le noyau de
//            Boussinesq sur les panneaux charges (methode de Steinbrenner),
//            collocation aux noeuds.
// Couplage : pression d'interaction constante par panneau tributaire Omega_j.
//
// Unites : kN, m  (=> pressions et modules en kPa).
//
// Traduction directe (numpy -> Eigen) du prototype Python Code_Python.cppnb.

#ifndef NOMINMAX
#define NOMINMAX  // evite que <windows.h> definisse des macros min/max qui entrent
                   // en conflit avec std::min/std::max et Eigen
#endif
#include <windows.h>
#include <shellapi.h>

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <tuple>
#include <utility>
#include <vector>

// ============================================================================
// 1. ELEMENT DE PLAQUE ACM
// ============================================================================
namespace {

    const Eigen::Vector4d XI(-1.0, 1.0, 1.0, -1.0);   // coins locaux 1..4, sens trigonometrique
    const Eigen::Vector4d ETA(-1.0, -1.0, 1.0, 1.0);

    // base polynomiale incomplete a 12 termes, en (xi, eta) :
    //   1, xi, eta, xi^2, xi*eta, eta^2, xi^3, xi^2*eta, xi*eta^2, eta^3, xi^3*eta, xi*eta^3
    const std::array<std::pair<int, int>, 12> EXP_ = { {
        {0, 0}, {1, 0}, {0, 1}, {2, 0}, {1, 1}, {0, 2},
        {3, 0}, {2, 1}, {1, 2}, {0, 3}, {3, 1}, {1, 3}
    } };

    double pow_(double v, int n) {
        return n < 0 ? 0.0 : std::pow(v, n);
    }

    // Vecteur (12) des monomes derives d_xi fois en xi et d_eta fois en eta.
    Eigen::VectorXd basis(double xi, double eta, int d_xi = 0, int d_eta = 0) {
        Eigen::VectorXd out = Eigen::VectorXd::Zero(12);
        for (int k = 0; k < 12; ++k) {
            const int p = EXP_[k].first, q = EXP_[k].second;
            double cp = 1.0;
            for (int i = 0; i < d_xi; ++i) cp *= (p - i);
            double cq = 1.0;
            for (int i = 0; i < d_eta; ++i) cq *= (q - i);
            if (cp == 0.0 || cq == 0.0) {
                out[k] = 0.0;
            }
            else {
                out[k] = cp * cq * pow_(xi, p - d_xi) * pow_(eta, q - d_eta);
            }
        }
        return out;
    }

    // Applique Pinv a un vecteur de base : equivalent numpy de (v @ Pinv).
    Eigen::VectorXd applyPinv(const Eigen::VectorXd& v, const Eigen::MatrixXd& Pinv) {
        return Pinv.transpose() * v;
    }

    // Inverse de la matrice 12x12 des conditions nodales (d = P c).
    Eigen::MatrixXd acmPinv(double a, double b) {
        Eigen::MatrixXd P(12, 12);
        for (int i = 0; i < 4; ++i) {
            P.row(3 * i) = basis(XI[i], ETA[i]).transpose();
            P.row(3 * i + 1) = (basis(XI[i], ETA[i], 1, 0) / a).transpose();  // ddl w_x = w_xi / a
            P.row(3 * i + 2) = (basis(XI[i], ETA[i], 0, 1) / b).transpose(); // ddl w_y = w_eta / b
        }
        return P.inverse();
    }

    // Les 12 fonctions de forme, ordre (w1,wx1,wy1,...,wy4).
    Eigen::VectorXd acmN(double xi, double eta, double a, double b,
        const Eigen::MatrixXd* PinvIn = nullptr) {
        Eigen::MatrixXd local;
        const Eigen::MatrixXd& Pinv = PinvIn ? *PinvIn : (local = acmPinv(a, b));
        return applyPinv(basis(xi, eta), Pinv);
    }

    // B (3x12) : kappa = {-w_xx, -w_yy, -2 w_xy} = B d.
    Eigen::MatrixXd acmB(double xi, double eta, double a, double b,
        const Eigen::MatrixXd* PinvIn = nullptr) {
        Eigen::MatrixXd local;
        const Eigen::MatrixXd& Pinv = PinvIn ? *PinvIn : (local = acmPinv(a, b));
        Eigen::VectorXd d2x = basis(xi, eta, 2, 0) / (a * a);
        Eigen::VectorXd d2y = basis(xi, eta, 0, 2) / (b * b);
        Eigen::VectorXd dxy = basis(xi, eta, 1, 1) / (a * b);
        Eigen::MatrixXd B(3, 12);
        B.row(0) = -applyPinv(d2x, Pinv).transpose();
        B.row(1) = -applyPinv(d2y, Pinv).transpose();
        B.row(2) = -applyPinv(2.0 * dxy, Pinv).transpose();
        return B;
    }

    // Matrice de comportement flexionnel (3x3) et rigidite D0.
    std::pair<Eigen::Matrix3d, double> dMatrix(double E, double nu, double h) {
        const double D0 = E * h * h * h / (12.0 * (1.0 - nu * nu));
        Eigen::Matrix3d D;
        D << 1.0, nu, 0.0,
            nu, 1.0, 0.0,
            0.0, 0.0, (1.0 - nu) / 2.0;
        D *= D0;
        return { D, D0 };
    }

    // Points et poids de Gauss-Legendre a n points sur [-1, 1] (Newton-Raphson).
    std::pair<Eigen::VectorXd, Eigen::VectorXd> gauss(int n) {
        Eigen::VectorXd x(n), w(n);
        const double PI = 3.14159265358979323846;
        for (int i = 0; i < n; ++i) {
            double xi = std::cos(PI * (i + 0.75) / (n + 0.5));
            double p1 = 0.0, dpn = 0.0;
            for (int iter = 0; iter < 100; ++iter) {
                double p0 = 1.0;
                p1 = xi;
                for (int k = 2; k <= n; ++k) {
                    double p2 = ((2 * k - 1) * xi * p1 - (k - 1) * p0) / k;
                    p0 = p1;
                    p1 = p2;
                }
                dpn = n * (xi * p1 - p0) / (xi * xi - 1.0);
                double dx = p1 / dpn;
                xi -= dx;
                if (std::abs(dx) < 1e-14) break;
            }
            x[i] = xi;
            w[i] = 2.0 / ((1.0 - xi * xi) * dpn * dpn);
        }
        std::vector<int> idx(n);
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b) { return x[a] < x[b]; });
        Eigen::VectorXd xs(n), ws(n);
        for (int i = 0; i < n; ++i) {
            xs[i] = x[idx[i]];
            ws[i] = w[idx[i]];
        }
        return { xs, ws };
    }

    // Matrice de rigidite elementaire 12x12 (integration de Gauss exacte).
    Eigen::MatrixXd keElement(double a, double b, double E, double nu, double h, int ngp = 4) {
        auto [D, D0] = dMatrix(E, nu, h);
        (void)D0;
        auto [gp, gw] = gauss(ngp);
        Eigen::MatrixXd Pinv = acmPinv(a, b);
        Eigen::MatrixXd K = Eigen::MatrixXd::Zero(12, 12);
        for (int i = 0; i < gp.size(); ++i) {
            for (int j = 0; j < gp.size(); ++j) {
                Eigen::MatrixXd B = acmB(gp[i], gp[j], a, b, &Pinv);
                K += gw[i] * gw[j] * (B.transpose() * D * B) * a * b;
            }
        }
        return K;
    }

    // Table canonique T[k, m] = integrale de la fonction locale k sur le quart m.
    Eigen::MatrixXd quarterTable(double a, double b, int ngp = 4) {
        auto [gp, gw] = gauss(ngp);
        Eigen::MatrixXd Pinv = acmPinv(a, b);
        Eigen::MatrixXd T = Eigen::MatrixXd::Zero(12, 4);
        // quart m : demi-domaine autour du coin m, en (xi,eta)
        const double bounds[4][4] = {
            {-1, 0, -1, 0},  // m=0
            {0, 1, -1, 0},   // m=1
            {0, 1, 0, 1},    // m=2
            {-1, 0, 0, 1}    // m=3
        };
        for (int m = 0; m < 4; ++m) {
            const double x1 = bounds[m][0], x2 = bounds[m][1];
            const double y1 = bounds[m][2], y2 = bounds[m][3];
            const double jx = (x2 - x1) / 2.0, cx = (x2 + x1) / 2.0;
            const double jy = (y2 - y1) / 2.0, cy = (y2 + y1) / 2.0;
            for (int i = 0; i < gp.size(); ++i) {
                for (int j = 0; j < gp.size(); ++j) {
                    T.col(m) += gw[i] * gw[j] *
                        acmN(cx + jx * gp[i], cy + jy * gp[j], a, b, &Pinv) *
                        jx * jy * a * b;
                }
            }
        }
        return T;
    }

}  // namespace

// Test d'appartenance point-polygone par lancer de rayon (ray casting) : compte le
// nombre d'intersections entre une demi-droite horizontale partant de (x,y) et les
// aretes du contour ; un nombre impair signifie que le point est a l'interieur.
// Fonctionne pour tout polygone simple (convexe ou non), quel que soit le sens de
// parcours des sommets -- c'est la brique de base qui permet au maillage de gerer
// aussi bien un rectangle qu'un radier en L, en U, ou toute autre forme.
bool pointDansPolygone(double x, double y, const std::vector<Eigen::Vector2d>& contour) {
    bool dedans = false;
    const int n = static_cast<int>(contour.size());
    for (int i = 0, j = n - 1; i < n; j = i++) {
        const double xi = contour[i].x(), yi = contour[i].y();
        const double xj = contour[j].x(), yj = contour[j].y();
        const bool traverse = (yi > y) != (yj > y);
        if (traverse && x < (xj - xi) * (y - yi) / (yj - yi) + xi) dedans = !dedans;
    }
    return dedans;
}

// Distance du point (px,py) au segment [P1,P2].
double distancePointSegment(double px, double py, double x1, double y1, double x2, double y2) {
    const double dx = x2 - x1, dy = y2 - y1;
    const double long2 = dx * dx + dy * dy;
    if (long2 < 1e-20) return std::hypot(px - x1, py - y1);
    double t = ((px - x1) * dx + (py - y1) * dy) / long2;
    t = std::clamp(t, 0.0, 1.0);
    return std::hypot(px - (x1 + t * dx), py - (y1 + t * dy));
}

// Comme pointDansPolygone, mais accepte aussi un point situe exactement sur le contour
// (bord ou sommet, y compris un coin du radier). Un lancer de rayon strict classe un
// point sur la frontiere de maniere ambigue (souvent "dehors"), ce qui rejetait a tort
// une charge placee volontairement en bordure ou dans un coin -- cas pourtant tres
// courant (charge de rive, poteau d'angle). Tolerance de 1 mm : suffisante pour couvrir
// les arrondis de saisie/flottants, sans accepter un point reellement hors du radier.
bool pointDansOuSurPolygone(double x, double y, const std::vector<Eigen::Vector2d>& contour) {
    if (pointDansPolygone(x, y, contour)) return true;
    constexpr double TOLERANCE_BORD_M = 1e-3;
    const int n = static_cast<int>(contour.size());
    for (int i = 0, j = n - 1; i < n; j = i++) {
        if (distancePointSegment(x, y, contour[j].x(), contour[j].y(), contour[i].x(),
            contour[i].y()) <= TOLERANCE_BORD_M) {
            return true;
        }
    }
    return false;
}

// ============================================================================
// 2. MAILLAGE (rectangle, ou toute forme polygonale : L, U, ...)
// ============================================================================
// Grille de fond reguliere nx x ny sur la boite englobante [0,Lx] x [0,Ly], dont on
// ne conserve que les elements dont le centre tombe a l'interieur du contour fourni
// (algorithme du rayon). Noeuds et elements actifs sont ensuite renumerotes de
// maniere compacte (sans trou), dans le meme ordre de balayage que l'ancien
// maillage plein (ligne par ligne, gauche a droite) : pour un contour rectangulaire
// couvrant exactement [0,Lx] x [0,Ly], toutes les cellules sont actives et la
// numerotation obtenue est rigoureusement identique a l'ancien maillage dense.
class Mesh {
public:
    double Lx, Ly, a, b;
    int nx, ny, nn, ndof;
    std::vector<Eigen::Vector2d> nodes;
    std::vector<std::array<int, 4>> elements;

    // grilleVersNoeud[j*(nx+1)+i] / grilleVersElement[j*nx+i] : indice compact du
    // noeud/element de la grille de fond en position (i,j), ou -1 s'il est hors du
    // domaine actif (ex. dans l'echancrure d'un L ou d'un U). elementIJ[e] retrouve
    // la position (i,j) d'un element actif e, utile pour la detection de bord.
    std::vector<int> grilleVersNoeud, grilleVersElement;
    std::vector<std::array<int, 2>> elementIJ;

    Mesh(double Lx_, double Ly_, int nx_, int ny_, const std::vector<Eigen::Vector2d>& contour)
        : Lx(Lx_), Ly(Ly_), nx(nx_), ny(ny_) {
        a = Lx / nx / 2.0;
        b = Ly / ny / 2.0;
        Eigen::VectorXd xs = Eigen::VectorXd::LinSpaced(nx + 1, 0.0, Lx);
        Eigen::VectorXd ys = Eigen::VectorXd::LinSpaced(ny + 1, 0.0, Ly);

        std::vector<char> celluleActive(static_cast<std::size_t>(nx) * ny, 0);
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const double xc = 0.5 * (xs[i] + xs[i + 1]), yc = 0.5 * (ys[j] + ys[j + 1]);
                celluleActive[j * nx + i] = pointDansPolygone(xc, yc, contour) ? 1 : 0;
            }
        }

        std::vector<char> noeudActif(static_cast<std::size_t>(nx + 1) * (ny + 1), 0);
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                if (!celluleActive[j * nx + i]) continue;
                noeudActif[j * (nx + 1) + i] = 1;
                noeudActif[j * (nx + 1) + i + 1] = 1;
                noeudActif[(j + 1) * (nx + 1) + i] = 1;
                noeudActif[(j + 1) * (nx + 1) + i + 1] = 1;
            }
        }

        grilleVersNoeud.assign(static_cast<std::size_t>(nx + 1) * (ny + 1), -1);
        for (int j = 0; j <= ny; ++j) {
            for (int i = 0; i <= nx; ++i) {
                const int idx = j * (nx + 1) + i;
                if (!noeudActif[idx]) continue;
                grilleVersNoeud[idx] = static_cast<int>(nodes.size());
                nodes.emplace_back(xs[i], ys[j]);
            }
        }
        nn = static_cast<int>(nodes.size());
        ndof = 3 * nn;

        grilleVersElement.assign(static_cast<std::size_t>(nx) * ny, -1);
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const int idxCell = j * nx + i;
                if (!celluleActive[idxCell]) continue;
                const int n1 = j * (nx + 1) + i;
                grilleVersElement[idxCell] = static_cast<int>(elements.size());
                elements.push_back({ grilleVersNoeud[n1], grilleVersNoeud[n1 + 1],
                                      grilleVersNoeud[n1 + nx + 2], grilleVersNoeud[n1 + nx + 1] });
                elementIJ.push_back({ i, j });
            }
        }
    }

    // Nombre total de cellules de la grille de fond (actives + eliminees).
    int nbCellulesGrille() const { return nx * ny; }
    // Nombre de cellules eliminees car hors du contour (echancrure du L/U).
    int nbCellulesEliminees() const {
        return nbCellulesGrille() - static_cast<int>(elements.size());
    }

    // La cellule de grille (i,j) existe-t-elle et est-elle active ?
    bool estActive(int i, int j) const {
        if (i < 0 || i >= nx || j < 0 || j >= ny) return false;
        return grilleVersElement[j * nx + i] >= 0;
    }

    // Les 12 ddl globaux de l'element e, dans l'ordre local.
    std::array<int, 12> dofs(int e) const {
        std::array<int, 12> d{};
        for (int k = 0; k < 4; ++k) {
            const int n = elements[e][k];
            d[3 * k] = 3 * n;
            d[3 * k + 1] = 3 * n + 1;
            d[3 * k + 2] = 3 * n + 2;
        }
        return d;
    }

    // Aire A_j du panneau tributaire de chaque noeud.
    Eigen::VectorXd panelArea() const {
        Eigen::VectorXd A = Eigen::VectorXd::Zero(nn);
        for (const auto& el : elements) {
            for (int n : el) A[n] += a * b;
        }
        return A;
    }

    // Rectangle enveloppe du panneau j : (x1, x2, y1, y2), borne a la boite
    // englobante [0,Lx] x [0,Ly]. Approximation : pour un noeud situe au coin d'une
    // echancrure (L/U), le panneau n'est pas decoupe exactement le long du contour
    // reel -- simplification usuelle et documentee, l'erreur ne portant que sur les
    // quelques noeuds adjacents a l'echancrure.
    std::array<double, 4> panelBox(int j) const {
        const double x = nodes[j].x(), y = nodes[j].y();
        return { std::max(0.0, x - a), std::min(Lx, x + a),
                std::max(0.0, y - b), std::min(Ly, y + b) };
    }

    // ---- Contours pratiques -------------------------------------------------
    static std::vector<Eigen::Vector2d> contourRectangle(double Lx_, double Ly_) {
        return { {0.0, 0.0}, {Lx_, 0.0}, {Lx_, Ly_}, {0.0, Ly_} };
    }

    // Radier en L : rectangle Lx x Ly auquel on retire un rectangle Bx x By dans le
    // coin superieur droit (0 < Bx < Lx, 0 < By < Ly).
    static std::vector<Eigen::Vector2d> contourL(double Lx_, double Ly_, double Bx, double By) {
        return { {0.0, 0.0}, {Lx_, 0.0}, {Lx_, Ly_ - By}, {Lx_ - Bx, Ly_ - By},
                 {Lx_ - Bx, Ly_}, {0.0, Ly_} };
    }

    // Radier en U : rectangle Lx x Ly auquel on retire une echancrure rectangulaire
    // en partie superieure, entre les abscisses x1 et x2 (0 <= x1 < x2 <= Lx), sur
    // une profondeur By depuis le haut (0 < By < Ly).
    static std::vector<Eigen::Vector2d> contourU(double Lx_, double Ly_, double x1, double x2,
        double By) {
        return { {0.0, 0.0}, {Lx_, 0.0}, {Lx_, Ly_}, {x2, Ly_}, {x2, Ly_ - By},
                 {x1, Ly_ - By}, {x1, Ly_}, {0.0, Ly_} };
    }
};

// ============================================================================
// 3. ASSEMBLAGE PLAQUE ET MATRICE DE PASSAGE
// ============================================================================
// Charge ponctuelle (kN) au point de coordonnees globales (x, y) (m) -- distribuee de maniere
// coherente (via les fonctions de forme) sur l'element qui contient ce point, plutot que
// d'exiger que le point tombe exactement sur un noeud.
struct PointLoad {
    double x, y, P;
};

// Charge surfacique (kPa) uniforme sur le rectangle [xmin,xmax] x [ymin,ymax] (m), qui peut
// ne recouvrir qu'une partie du radier.
struct SurfaceLoad {
    double xmin, xmax, ymin, ymax, q;
};

// Charge lineique (kN/m) le long du segment (x1,y1)-(x2,y2), en coordonnees globales (m).
struct LineLoad {
    double x1, y1, x2, y2, q;
};

// Clip de Liang-Barsky : intersecte le segment (x1,y1)-(x2,y2) avec le rectangle
// [xmin,xmax] x [ymin,ymax]. Si intersection non vide, renvoie true et les parametres
// t0 <= t1 dans [0,1] du sous-segment inclus dans le rectangle.
//
// inclureBordXMax/inclureBordYMax : quand un segment est parallele a un bord et tombe
// exactement sur une frontiere partagee par deux elements (ex. une charge lineique
// horizontale qui coincide avec une ligne du maillage), il ne doit etre compte que par
// l'un des deux elements adjacents, sinon la charge equivalente est comptee deux fois. Par
// convention, chaque element accepte son bord xmin/ymin, mais pas son bord xmax/ymax --
// sauf s'il s'agit du bord exterieur du radier (dernier element de la rangee/colonne), qui
// doit rester inclusif pour ne pas perdre la charge quand elle coincide avec ce bord.
bool clipSegmentToBox(double x1, double y1, double x2, double y2, double xmin, double xmax,
    double ymin, double ymax, bool inclureBordXMax, bool inclureBordYMax,
    double& t0, double& t1) {
    const double dx = x2 - x1, dy = y2 - y1;
    double tE = 0.0, tL = 1.0;
    auto clipTest = [&](double p, double q, bool inclusifSiParallele) {
        if (std::abs(p) < 1e-14) return inclusifSiParallele ? (q >= 0.0) : (q > 0.0);
        const double r = q / p;
        if (p < 0.0) {
            if (r > tL) return false;
            if (r > tE) tE = r;
        }
        else {
            if (r < tE) return false;
            if (r < tL) tL = r;
        }
        return true;
        };
    if (!clipTest(-dx, x1 - xmin, true)) return false;
    if (!clipTest(dx, xmax - x1, inclureBordXMax)) return false;
    if (!clipTest(-dy, y1 - ymin, true)) return false;
    if (!clipTest(dy, ymax - y1, inclureBordYMax)) return false;
    if (tE > tL) return false;
    t0 = tE;
    t1 = tL;
    return true;
}

// Vecteur de charge equivalent global pour des charges lineiques : chaque segment est
// decoupe (clip) element par element, puis integre par quadrature de Gauss 1D le long du
// sous-segment inclus dans l'element, avec les fonctions de forme de la plaque evaluees en
// chaque point d'integration (charge repartie de maniere coherente, comme feElement pour une
// pression surfacique).
Eigen::VectorXd assembleFLignes(const Mesh& mesh, const std::vector<LineLoad>& lineLoads,
    int ngp = 4) {
    Eigen::VectorXd F = Eigen::VectorXd::Zero(mesh.ndof);
    if (lineLoads.empty()) return F;
    auto [gp, gw] = gauss(ngp);
    const Eigen::MatrixXd Pinv = acmPinv(mesh.a, mesh.b);  // uniforme : maillage regulier

    for (const auto& ll : lineLoads) {
        for (std::size_t e = 0; e < mesh.elements.size(); ++e) {
            const auto& el = mesh.elements[e];
            const Eigen::Vector2d& c0 = mesh.nodes[el[0]];  // coin inferieur gauche
            const Eigen::Vector2d& c2 = mesh.nodes[el[2]];  // coin superieur droit
            const double xmin = c0.x(), ymin = c0.y(), xmax = c2.x(), ymax = c2.y();
            // Un bord est "exterieur" (traitement inclusif, cf. clipSegmentToBox) des qu'il
            // n'a pas de voisin actif de ce cote -- que ce soit le bord de la boite englobante
            // ou le bord d'une echancrure interne (radier en L/U).
            const auto [gi, gj] = mesh.elementIJ[e];
            const bool bordDroit = !mesh.estActive(gi + 1, gj);
            const bool bordHaut = !mesh.estActive(gi, gj + 1);

            double t0, t1;
            if (!clipSegmentToBox(ll.x1, ll.y1, ll.x2, ll.y2, xmin, xmax, ymin, ymax, bordDroit,
                bordHaut, t0, t1))
                continue;

            const double xA = ll.x1 + t0 * (ll.x2 - ll.x1), yA = ll.y1 + t0 * (ll.y2 - ll.y1);
            const double xB = ll.x1 + t1 * (ll.x2 - ll.x1), yB = ll.y1 + t1 * (ll.y2 - ll.y1);
            const double demiLongueur = 0.5 * std::sqrt((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA));
            if (demiLongueur < 1e-14) continue;  // le segment ne fait qu'effleurer l'element

            const double xc = 0.5 * (xmin + xmax), yc = 0.5 * (ymin + ymax);
            const auto d = mesh.dofs(static_cast<int>(e));

            for (int k = 0; k < gp.size(); ++k) {
                const double s = gp[k];
                const double x = 0.5 * (xA + xB) + 0.5 * s * (xB - xA);
                const double y = 0.5 * (yA + yB) + 0.5 * s * (yB - yA);
                const double xi = (x - xc) / mesh.a;
                const double eta = (y - yc) / mesh.b;
                const Eigen::VectorXd N = acmN(xi, eta, mesh.a, mesh.b, &Pinv);
                for (int i = 0; i < 12; ++i) F[d[i]] += gw[k] * ll.q * N[i] * demiLongueur;
            }
        }
    }
    return F;
}

Eigen::MatrixXd assembleK(const Mesh& mesh, double E, double nu, double h) {
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(mesh.ndof, mesh.ndof);
    Eigen::MatrixXd ke = keElement(mesh.a, mesh.b, E, nu, h);
    for (std::size_t e = 0; e < mesh.elements.size(); ++e) {
        auto d = mesh.dofs(static_cast<int>(e));
        for (int i = 0; i < 12; ++i)
            for (int j = 0; j < 12; ++j)
                K(d[i], d[j]) += ke(i, j);
    }
    return K;
}

// Localise l'element du maillage qui contient (xp, yp) (bornes de la boite englobante
// incluses) et renvoie ses coordonnees locales (xi, eta) dans [-1, 1]. Si (xp, yp) tombe
// dans une zone hors du domaine actif (ex. l'echancrure d'un radier en L/U), la charge
// est rattachee a l'element actif le plus proche plutot que d'etre perdue ou de planter.
void localiserPoint(const Mesh& mesh, double xp, double yp, int& e, double& xi, double& eta) {
    int i = std::clamp(static_cast<int>(xp / (mesh.Lx / mesh.nx)), 0, mesh.nx - 1);
    int j = std::clamp(static_cast<int>(yp / (mesh.Ly / mesh.ny)), 0, mesh.ny - 1);

    if (mesh.grilleVersElement[j * mesh.nx + i] < 0) {
        int meilleurI = -1, meilleurJ = -1;
        double meilleureDist = std::numeric_limits<double>::infinity();
        for (int jj = 0; jj < mesh.ny; ++jj) {
            for (int ii = 0; ii < mesh.nx; ++ii) {
                if (mesh.grilleVersElement[jj * mesh.nx + ii] < 0) continue;
                const double dx = static_cast<double>(ii - i), dy = static_cast<double>(jj - j);
                const double d = dx * dx + dy * dy;
                if (d < meilleureDist) { meilleureDist = d; meilleurI = ii; meilleurJ = jj; }
            }
        }
        std::cerr << "Attention : le point (" << xp << ", " << yp << ") tombe hors du domaine "
            "actif du radier (echancrure) ; la charge est rattachee a l'element actif le plus "
            "proche.\n";
        i = meilleurI;
        j = meilleurJ;
    }

    e = mesh.grilleVersElement[j * mesh.nx + i];
    const auto& el = mesh.elements[e];
    const Eigen::Vector2d& c0 = mesh.nodes[el[0]];
    const Eigen::Vector2d& c2 = mesh.nodes[el[2]];
    const double xc = 0.5 * (c0.x() + c2.x()), yc = 0.5 * (c0.y() + c2.y());
    xi = (xp - xc) / mesh.a;
    eta = (yp - yc) / mesh.b;
}

// Vecteur de charge equivalent global pour des charges ponctuelles definies par leurs
// coordonnees : chaque force est distribuee de maniere coherente (fonctions de forme de la
// plaque) sur l'element qui contient son point d'application.
Eigen::VectorXd assembleFPoints(const Mesh& mesh, const std::vector<PointLoad>& pointLoads) {
    Eigen::VectorXd F = Eigen::VectorXd::Zero(mesh.ndof);
    if (pointLoads.empty()) return F;
    const Eigen::MatrixXd Pinv = acmPinv(mesh.a, mesh.b);
    for (const auto& pl : pointLoads) {
        int e;
        double xi, eta;
        localiserPoint(mesh, pl.x, pl.y, e, xi, eta);
        const Eigen::VectorXd N = acmN(xi, eta, mesh.a, mesh.b, &Pinv);
        const auto d = mesh.dofs(e);
        for (int k = 0; k < 12; ++k) F[d[k]] += pl.P * N[k];
    }
    return F;
}

// Vecteur de charge equivalent global pour des charges surfaciques bornees : chaque
// rectangle charge est decoupe (intersecte) element par element, puis integre par
// quadrature de Gauss sur la portion effectivement recouverte (charge repartie de maniere
// coherente, comme pour une pression uniforme sur tout un element, mais restreinte a la
// zone definie par l'utilisateur).
Eigen::VectorXd assembleFSurfaces(const Mesh& mesh, const std::vector<SurfaceLoad>& surfaceLoads,
    int ngp = 4) {
    Eigen::VectorXd F = Eigen::VectorXd::Zero(mesh.ndof);
    if (surfaceLoads.empty()) return F;
    auto [gp, gw] = gauss(ngp);
    const Eigen::MatrixXd Pinv = acmPinv(mesh.a, mesh.b);  // uniforme : maillage regulier

    for (const auto& sl : surfaceLoads) {
        for (std::size_t e = 0; e < mesh.elements.size(); ++e) {
            const auto& el = mesh.elements[e];
            const Eigen::Vector2d& c0 = mesh.nodes[el[0]];
            const Eigen::Vector2d& c2 = mesh.nodes[el[2]];
            const double ex0 = c0.x(), ey0 = c0.y(), ex1 = c2.x(), ey1 = c2.y();

            const double ox0 = std::max(ex0, sl.xmin), ox1 = std::min(ex1, sl.xmax);
            const double oy0 = std::max(ey0, sl.ymin), oy1 = std::min(ey1, sl.ymax);
            if (ox1 <= ox0 || oy1 <= oy0) continue;  // pas de recouvrement avec cet element

            const double xc = 0.5 * (ex0 + ex1), yc = 0.5 * (ey0 + ey1);
            const double xi0 = (ox0 - xc) / mesh.a, xi1 = (ox1 - xc) / mesh.a;
            const double eta0 = (oy0 - yc) / mesh.b, eta1 = (oy1 - yc) / mesh.b;
            const double jxi = 0.5 * (xi1 - xi0), cxi = 0.5 * (xi1 + xi0);
            const double jeta = 0.5 * (eta1 - eta0), ceta = 0.5 * (eta1 + eta0);

            const auto d = mesh.dofs(static_cast<int>(e));
            for (int i = 0; i < gp.size(); ++i) {
                for (int j = 0; j < gp.size(); ++j) {
                    const Eigen::VectorXd N =
                        acmN(cxi + jxi * gp[i], ceta + jeta * gp[j], mesh.a, mesh.b, &Pinv);
                    const double poids = gw[i] * gw[j] * sl.q * jxi * jeta * mesh.a * mesh.b;
                    for (int k = 0; k < 12; ++k) F[d[k]] += poids * N[k];
                }
            }
        }
    }
    return F;
}

// Mp[i, j] = integrale sur Omega_j de N_i dA   (3n x n).
Eigen::MatrixXd assembleMp(const Mesh& mesh) {
    Eigen::MatrixXd T = quarterTable(mesh.a, mesh.b);
    Eigen::MatrixXd Mp = Eigen::MatrixXd::Zero(mesh.ndof, mesh.nn);
    for (std::size_t e = 0; e < mesh.elements.size(); ++e) {
        auto d = mesh.dofs(static_cast<int>(e));
        for (int m = 0; m < 4; ++m) {              // quart du coin local m
            const int nj = mesh.elements[e][m];
            for (int row = 0; row < 12; ++row) Mp(d[row], nj) += T(row, m);
        }
    }
    return Mp;
}

// ============================================================================
// 4. SOL : NOYAU DE BOUSSINESQ INTEGRE (STEINBRENNER)
// ============================================================================
// Tassement vertical a la profondeur z sous l'angle d'un rectangle A x B
// charge par une pression unitaire, demi-espace elastique (E, nu).
double sCorner(double A, double B, double z, double E, double nu) {
    A = std::abs(A);
    B = std::abs(B);
    if (A < 1e-14 || B < 1e-14) return 0.0;
    const double R = std::sqrt(A * A + B * B + z * z);
    double J1, J2;
    if (z < 1e-12) {
        J1 = A * std::log((B + std::sqrt(A * A + B * B)) / A) +
            B * std::log((A + std::sqrt(A * A + B * B)) / B);
        J2 = 0.0;
    }
    else {
        J1 = A * std::log((B + R) / std::sqrt(A * A + z * z)) +
            B * std::log((A + R) / std::sqrt(B * B + z * z)) -
            z * std::atan(A * B / (z * R));
        J2 = z * std::atan(A * B / (z * R));
    }
    const double PI = 3.14159265358979323846;
    return (1.0 + nu) / (2.0 * PI * E) * (J2 + 2.0 * (1.0 - nu) * J1);
}

double sign_(double v) {
    return (v > 0.0) - (v < 0.0);
}

// Extension impaire du terme d'angle : rectangles signes.
double phi_(double u, double v, double z, double E, double nu) {
    return sign_(u) * sign_(v) * sCorner(u, v, z, E, nu);
}

struct Layer {
    double zt, zb, E, nu;  // zb = infini pour un demi-espace
};

// Tassement en surface au point (xp, yp) du a une pression unitaire
// uniforme sur le rectangle box = (x1, x2, y1, y2).
double settlement(double xp, double yp, const std::array<double, 4>& box,
    const std::vector<Layer>& layers) {
    const double x1 = box[0], x2 = box[1], y1 = box[2], y2 = box[3];
    const double u1 = x1 - xp, u2 = x2 - xp;
    const double v1 = y1 - yp, v2 = y2 - yp;
    double s = 0.0;
    for (const auto& layer : layers) {
        for (auto [zz, sg] : { std::pair{layer.zt, 1.0}, std::pair{layer.zb, -1.0} }) {
            if (!std::isfinite(zz)) continue;
            s += sg * (phi_(u2, v2, zz, layer.E, layer.nu) - phi_(u1, v2, zz, layer.E, layer.nu) -
                phi_(u2, v1, zz, layer.E, layer.nu) + phi_(u1, v1, zz, layer.E, layer.nu));
        }
    }
    return s;
}

// Ts[i, j] = tassement au noeud i sous pression unitaire sur Omega_j.
Eigen::MatrixXd assembleTs(const Mesh& mesh, const std::vector<Layer>& layers) {
    Eigen::MatrixXd Ts = Eigen::MatrixXd::Zero(mesh.nn, mesh.nn);
    for (int j = 0; j < mesh.nn; ++j) {
        auto box = mesh.panelBox(j);
        for (int i = 0; i < mesh.nn; ++i) {
            Ts(i, j) = settlement(mesh.nodes[i].x(), mesh.nodes[i].y(), box, layers);
        }
    }
    return Ts;
}

// ============================================================================
// 5. SYSTEME COUPLE
// ============================================================================
struct PlateProps {
    double E, nu, h;
};

struct SsiResult {
    Eigen::MatrixXd K, Mp, Ts, Shat, Ks, G, Ksol, Kg;
    Eigen::VectorXd F, A, a, sbar, P, p;
};

SsiResult solveSsi(const Mesh& mesh, const PlateProps& plate, const std::vector<Layer>& layers,
    const std::vector<SurfaceLoad>& surfaceLoads = {},
    const std::vector<PointLoad>& pointLoads = {},
    const std::vector<LineLoad>& lineLoads = {}, bool symmetrize = true) {
    SsiResult r;
    r.K = assembleK(mesh, plate.E, plate.nu, plate.h);
    r.F = assembleFSurfaces(mesh, surfaceLoads) + assembleFPoints(mesh, pointLoads) +
        assembleFLignes(mesh, lineLoads);
    r.Mp = assembleMp(mesh);
    r.Ts = assembleTs(mesh, layers);
    r.A = mesh.panelArea();

    r.Shat = r.Ts;
    for (int j = 0; j < mesh.nn; ++j) r.Shat.col(j) /= r.A[j];  // souplesse vis-a-vis des forces de panneau
    if (symmetrize) r.Shat = 0.5 * (r.Shat + r.Shat.transpose());
    r.Ks = r.Shat.inverse();  // rigidite du sol (forces de panneau)

    Eigen::MatrixXd MpOverA = r.Mp;
    for (int j = 0; j < mesh.nn; ++j) MpOverA.col(j) /= r.A[j];
    r.G = MpOverA.transpose();  // s_bar = G a

    r.Ksol = r.G.transpose() * r.Ks * r.G;
    r.Kg = r.K + r.Ksol;
    r.a = r.Kg.colPivHouseholderQr().solve(r.F);
    r.sbar = r.G * r.a;
    r.P = r.Ks * r.sbar;                 // forces de panneau
    r.p = r.P.array() / r.A.array();     // pressions d'interaction
    return r;
}

// Moments {Mx, My, Mxy} au point (xi, eta) de l'element e.
Eigen::Vector3d elementMoments(const Mesh& mesh, const Eigen::VectorXd& aSol,
    const PlateProps& plate, int e, double xi, double eta) {
    auto [D, D0] = dMatrix(plate.E, plate.nu, plate.h);
    (void)D0;
    Eigen::MatrixXd B = acmB(xi, eta, mesh.a, mesh.b);
    auto d = mesh.dofs(e);
    Eigen::VectorXd de(12);
    for (int i = 0; i < 12; ++i) de[i] = aSol[d[i]];
    return D * B * de;
}

// ============================================================================
// Saisie interactive et execution
// (le prototype Python original ne contenait que des definitions, sans point
// d'entree ; ce main() lit tous les parametres necessaires au calcul).
//
// Chaque invite affiche une valeur par defaut entre crochets : appuyer sur
// Entree sans rien taper conserve cette valeur, ce qui permet de relancer le
// calcul rapidement sans tout ressaisir a chaque fois.
// ============================================================================
namespace {

    // Lit une ligne complete sur stdin ; une ligne vide conserve la valeur par
    // defaut. On utilise getline (pas operator>>) pour que "Entree" seule soit
    // bien detectee comme choix par defaut plutot que de bloquer en attente
    // d'un prochain jeton.
    //
    // Si l'entree standard est fermee (EOF, ex. execution avec un flux redirige
    // qui s'epuise), on arrete proprement plutot que de renvoyer une chaine vide
    // indefiniment : sans ce garde-fou, une boucle de validation dont la valeur
    // par defaut ne satisfait jamais le critere (ex. un rectangle de charge nul
    // qui echoue toujours au test xmax > xmin) tournerait indefiniment.
    std::string lireLigne(const std::string& invite) {
        std::cout << invite;
        std::cout.flush();
        std::string ligne;
        if (!std::getline(std::cin, ligne)) {
            std::cerr << "\nEntree standard fermee (EOF) avant la fin de la saisie -- arret.\n";
            std::exit(1);
        }
        return ligne;
    }

    // Affiche une matrice complete, indices 1-based, notation scientifique.
    void afficherMatrice(const std::string& nom, const Eigen::MatrixXd& M) {
        std::cout << "\n----- " << nom << " (" << M.rows() << " x " << M.cols() << ") -----\n";
        std::cout << std::scientific << std::setprecision(4);
        std::cout << std::setw(8) << "";
        for (int j = 0; j < M.cols(); ++j) std::cout << std::setw(13) << (j + 1);
        std::cout << "\n";
        for (int i = 0; i < M.rows(); ++i) {
            std::cout << std::setw(8) << (i + 1);
            for (int j = 0; j < M.cols(); ++j) std::cout << std::setw(13) << M(i, j);
            std::cout << "\n";
        }
        std::cout << std::fixed << std::setprecision(6);
    }

    constexpr int LARGEUR_COL_NOEUD = 8;
    constexpr int LARGEUR_COL_VALEUR = 23;

    // Trace une ligne de separation "=" de la largeur totale du tableau (nbColonnes
    // colonnes de valeurs en plus de la colonne "Noeud").
    void separateurTableau(int nbColonnes) {
        std::cout << std::string(LARGEUR_COL_NOEUD + nbColonnes * LARGEUR_COL_VALEUR, '=') << "\n";
    }

    // Tableau noeud x 3 ddl (deplacements ou reactions du sol), notation
    // scientifique a 6 chiffres significatifs, colonnes alignees.
    void afficherTableauVecteur3Ddl(const std::string& titre,
        const std::array<std::string, 3>& colonnes, int nn, const Eigen::VectorXd& v) {
        std::cout << "\n===== " << titre << " (" << (3 * nn) << " x 1, " << nn << " noeuds) =====\n";
        separateurTableau(3);
        std::cout << std::setw(LARGEUR_COL_NOEUD) << "Noeud";
        for (const auto& c : colonnes) std::cout << std::setw(LARGEUR_COL_VALEUR) << c;
        std::cout << "\n";
        separateurTableau(3);
        std::cout << std::scientific << std::setprecision(6);
        for (int i = 0; i < nn; ++i) {
            std::cout << std::setw(LARGEUR_COL_NOEUD) << (i + 1);
            for (int k = 0; k < 3; ++k) std::cout << std::setw(LARGEUR_COL_VALEUR) << v[3 * i + k];
            std::cout << "\n";
        }
        separateurTableau(3);
        std::cout << std::fixed << std::setprecision(6);
    }

    // Tableau noeud x 1 valeur (pression ou coefficient de reaction), notation
    // scientifique a 6 chiffres significatifs, colonnes alignees.
    void afficherTableauVecteurNoeud(const std::string& titre, const std::string& colonne,
        int nn, const Eigen::VectorXd& v) {
        std::cout << "\n===== " << titre << " (" << nn << " x 1) =====\n";
        separateurTableau(1);
        std::cout << std::setw(LARGEUR_COL_NOEUD) << "Noeud" << std::setw(LARGEUR_COL_VALEUR) << colonne
            << "\n";
        separateurTableau(1);
        std::cout << std::scientific << std::setprecision(6);
        for (int i = 0; i < nn; ++i) {
            std::cout << std::setw(LARGEUR_COL_NOEUD) << (i + 1) << std::setw(LARGEUR_COL_VALEUR) << v[i]
                << "\n";
        }
        separateurTableau(1);
        std::cout << std::fixed << std::setprecision(6);
    }

    // ks(x,y) = p(x,y) / s(x,y), avec securite contre la division par un
    // tassement nul ou quasi nul (retourne 0 dans ce cas plutot que +-inf/NaN).
    Eigen::VectorXd coefficientReactionSol(const Eigen::VectorXd& p, const Eigen::VectorXd& s) {
        constexpr double SEUIL_TASSEMENT_NUL = 1e-12;
        Eigen::VectorXd ks = Eigen::VectorXd::Zero(p.size());
        for (int i = 0; i < p.size(); ++i) {
            if (std::abs(s[i]) > SEUIL_TASSEMENT_NUL) ks[i] = p[i] / s[i];
        }
        return ks;
    }

    // Palette de type "jet" (bleu -> cyan -> vert -> jaune -> rouge) pour t dans [0,1] : bleu =
    // valeur la plus faible, rouge = valeur la plus forte (ici, le tassement le plus marque).
    void couleurJet(double t, int& r, int& g, int& b) {
        t = std::max(0.0, std::min(1.0, t));
        double v = t * 4.0;
        double rf = std::max(0.0, std::min(1.0, std::min(v - 1.5, -v + 4.5)));
        double gf = std::max(0.0, std::min(1.0, std::min(v - 0.5, -v + 3.5)));
        double bf = std::max(0.0, std::min(1.0, std::min(v + 0.5, -v + 2.5)));
        r = static_cast<int>(std::round(rf * 255.0));
        g = static_cast<int>(std::round(gf * 255.0));
        b = static_cast<int>(std::round(bf * 255.0));
    }

    // Exporte une carte d'isovaleurs (SVG) du tassement du radier = deplacement vertical w de la
    // plaque (meme principe que exporterIsovaleursTassementSVG de DimRadierSol.cpp). Chaque element
    // est sous-decoupe en NSUB x NSUB petits quadrilateres, colores par interpolation bilineaire de
    // w a partir des 4 noeuds de coin. Degrade bleu (le moins tasse) -> rouge (le plus tasse), avec
    // legende graduee. Affiche le chemin du fichier genere a chaque appel, pour consultation.
    void exporterIsovaleursTassementSVG(const Mesh& mesh, const Eigen::VectorXd& aSol,
        const std::string& cheminFichier) {
        const double echelle = 100.0;
        const double margeGauche = 90.0;
        const double margeDroite = 140.0;  // place pour la legende verticale
        const double margeHaut = 60.0;
        const double margeBas = 90.0;

        const double largeurDessin = mesh.Lx * echelle;
        const double hauteurDessin = mesh.Ly * echelle;
        const double largeurTotale = largeurDessin + margeGauche + margeDroite;
        const double hauteurTotale = hauteurDessin + margeHaut + margeBas;

        auto toPxX = [&](double x) { return margeGauche + x * echelle; };
        auto toPxY = [&](double y) { return margeHaut + (mesh.Ly - y) * echelle; };

        double wMin = aSol[0], wMax = aSol[0];
        for (int n = 0; n < mesh.nn; ++n) {
            double w = aSol[3 * n];
            wMin = std::min(wMin, w);
            wMax = std::max(wMax, w);
        }
        double ecart = wMax - wMin;
        if (ecart < 1e-15) ecart = 1.0;  // evite une division par zero si le champ est uniforme

        std::ofstream f(cheminFichier);
        if (!f.is_open()) {
            std::cerr << "Erreur : impossible de creer le fichier " << cheminFichier << std::endl;
            return;
        }

        f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
        f << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << largeurTotale << "\" height=\""
            << hauteurTotale << "\" viewBox=\"0 0 " << largeurTotale << " " << hauteurTotale
            << "\" font-family=\"Arial, sans-serif\">\n";
        f << "  <rect x=\"0\" y=\"0\" width=\"" << largeurTotale << "\" height=\"" << hauteurTotale
            << "\" fill=\"white\"/>\n";
        f << "  <text x=\"" << margeGauche
            << "\" y=\"24\" font-size=\"15\" font-weight=\"bold\">Tassement du radier (w)</text>\n";

        // --- Carte de couleurs : sous-decoupage fin de chaque element, interpolation bilineaire ---
        const int NSUB = 12;
        f << "  <g stroke=\"none\">\n";
        for (const auto& el : mesh.elements) {
            const Eigen::Vector2d& n0 = mesh.nodes[el[0]];
            const Eigen::Vector2d& n1 = mesh.nodes[el[1]];
            const Eigen::Vector2d& n2 = mesh.nodes[el[2]];
            const Eigen::Vector2d& n3 = mesh.nodes[el[3]];
            const double w0 = aSol[3 * el[0]], w1 = aSol[3 * el[1]];
            const double w2 = aSol[3 * el[2]], w3 = aSol[3 * el[3]];

            auto position = [&](double xi, double eta) {
                double x = (1 - xi) * (1 - eta) * n0.x() + xi * (1 - eta) * n1.x() +
                    xi * eta * n2.x() + (1 - xi) * eta * n3.x();
                double y = (1 - xi) * (1 - eta) * n0.y() + xi * (1 - eta) * n1.y() +
                    xi * eta * n2.y() + (1 - xi) * eta * n3.y();
                return std::make_pair(x, y);
                };
            auto valeurW = [&](double xi, double eta) {
                return (1 - xi) * (1 - eta) * w0 + xi * (1 - eta) * w1 + xi * eta * w2 +
                    (1 - xi) * eta * w3;
                };

            for (int i = 0; i < NSUB; ++i) {
                for (int j = 0; j < NSUB; ++j) {
                    const double xi0 = static_cast<double>(i) / NSUB, xi1 = static_cast<double>(i + 1) / NSUB;
                    const double eta0 = static_cast<double>(j) / NSUB, eta1 = static_cast<double>(j + 1) / NSUB;

                    auto pa = position(xi0, eta0);
                    auto pb = position(xi1, eta0);
                    auto pc = position(xi1, eta1);
                    auto pd = position(xi0, eta1);

                    const double wCentre = valeurW((xi0 + xi1) / 2.0, (eta0 + eta1) / 2.0);
                    // Le tassement est negatif (vers le bas) : la valeur la plus negative (wMin) est
                    // le tassement le plus severe et doit ressortir en rouge (t=1), tandis que la
                    // valeur la plus grande (wMax, la moins tassee) doit ressortir en bleu (t=0).
                    const double t = (wMax - wCentre) / ecart;
                    int r, g, b;
                    couleurJet(t, r, g, b);

                    f << "    <polygon points=\"" << toPxX(pa.first) << "," << toPxY(pa.second) << " "
                        << toPxX(pb.first) << "," << toPxY(pb.second) << " " << toPxX(pc.first) << ","
                        << toPxY(pc.second) << " " << toPxX(pd.first) << "," << toPxY(pd.second)
                        << "\" fill=\"rgb(" << r << "," << g << "," << b << ")\"/>\n";
                }
            }
        }
        f << "  </g>\n";

        // --- Contour des elements (par-dessus la carte de couleurs) ---
        f << "  <g stroke=\"#222222\" stroke-width=\"1.2\" fill=\"none\">\n";
        for (const auto& el : mesh.elements) {
            f << "    <polygon points=\"";
            for (int k = 0; k < 4; ++k) {
                const Eigen::Vector2d& n = mesh.nodes[el[k]];
                f << toPxX(n.x()) << "," << toPxY(n.y());
                if (k < 3) f << " ";
            }
            f << "\"/>\n";
        }
        f << "  </g>\n";

        // --- Cotation X/Y ---
        const double yCoteX = margeHaut + hauteurDessin + 30;
        f << "  <text x=\"" << (margeGauche + largeurDessin / 2) << "\" y=\"" << (yCoteX + 18)
            << "\" text-anchor=\"middle\" font-size=\"11\">Lx = " << mesh.Lx << " m (axe X)</text>\n";

        const double xCoteY = margeGauche - 45;
        f << "  <text x=\"" << xCoteY << "\" y=\"" << (margeHaut + hauteurDessin / 2)
            << "\" text-anchor=\"middle\" font-size=\"11\" transform=\"rotate(-90 " << xCoteY << ","
            << (margeHaut + hauteurDessin / 2) << ")\">Ly = " << mesh.Ly << " m (axe Y)</text>\n";

        // --- Legende : degrade vertical bleu (bas) -> rouge (haut) avec graduations de w ---
        const double legX = margeGauche + largeurDessin + 30;
        const double legYHaut = margeHaut;
        const double legHauteur = hauteurDessin;
        const double legLargeur = 22.0;

        f << "  <text x=\"" << legX << "\" y=\"" << (legYHaut - 8) << "\" font-size=\"10\">Tassement w (m)</text>\n";

        const int NBANDES = 20;
        for (int k = 0; k < NBANDES; ++k) {
            const double tHaut = 1.0 - static_cast<double>(k) / NBANDES;
            const double tBas = 1.0 - static_cast<double>(k + 1) / NBANDES;
            const double tMil = (tHaut + tBas) / 2.0;
            int r, g, b;
            couleurJet(tMil, r, g, b);
            const double yHautPx = legYHaut + legHauteur * (static_cast<double>(k) / NBANDES);
            const double yBasPx = legYHaut + legHauteur * (static_cast<double>(k + 1) / NBANDES);
            f << "    <rect x=\"" << legX << "\" y=\"" << yHautPx << "\" width=\"" << legLargeur
                << "\" height=\"" << (yBasPx - yHautPx) << "\" fill=\"rgb(" << r << "," << g << "," << b
                << ")\"/>\n";

            // tHaut=1 (haut de legende, rouge) doit afficher wMin (tassement le plus severe) ;
            // tHaut=0 (bas de legende, bleu) doit afficher wMax (le moins tasse).
            const double wGrad = wMax - ecart * tHaut;
            f << "    <text x=\"" << (legX + legLargeur + 6) << "\" y=\"" << (yHautPx + 4)
                << "\" font-size=\"8.5\">" << std::scientific << std::setprecision(2) << wGrad
                << std::defaultfloat << "</text>\n";
        }
        f << "    <text x=\"" << (legX + legLargeur + 6) << "\" y=\"" << (legYHaut + legHauteur + 4)
            << "\" font-size=\"8.5\">" << std::scientific << std::setprecision(2) << wMax
            << std::defaultfloat << "</text>\n";
        f << "    <rect x=\"" << legX << "\" y=\"" << legYHaut << "\" width=\"" << legLargeur
            << "\" height=\"" << legHauteur << "\" fill=\"none\" stroke=\"black\" stroke-width=\"0.8\"/>\n";

        f << "</svg>\n";
        f.close();

        std::cout << "\nCarte d'isovaleurs du tassement exportee (format vectoriel SVG) : "
            << cheminFichier << std::endl;
    }

    // Exporte une page HTML autonome avec un graphique 3D interactif (Plotly.js, via CDN) du
    // tassement du radier : la surface z = tassement amplifie (pour rester visible a l'echelle du
    // radier), colorisee par la valeur reelle (rouge = tassement le plus severe, bleu = le plus
    // faible, meme convention que la carte SVG), avec le maillage EF trace sur la surface
    // (contours suivant exactement les lignes de noeuds). Affiche le chemin du fichier a chaque
    // appel, pour ouverture immediate dans un navigateur.
    void exporterGraphique3DTassementHTML(const Mesh& mesh, const Eigen::VectorXd& aSol,
        const std::string& cheminFichier) {
        // xs/ys sont recalcules directement sur la boite englobante (independants de l'ordre
        // de stockage, compacte, des noeuds actifs). Les positions de grille hors du domaine
        // actif (echancrure d'un radier en L/U) recoivent NaN, ecrit "null" en JS plus bas :
        // Plotly interprete cela comme un trou dans la surface plutot que comme un tassement
        // nul fictif.
        std::vector<double> xs(mesh.nx + 1), ys(mesh.ny + 1);
        for (int i = 0; i <= mesh.nx; ++i) xs[i] = i * mesh.Lx / mesh.nx;
        for (int j = 0; j <= mesh.ny; ++j) ys[j] = j * mesh.Ly / mesh.ny;

        const double nan = std::numeric_limits<double>::quiet_NaN();
        std::vector<std::vector<double>> Z(mesh.ny + 1, std::vector<double>(mesh.nx + 1, nan));
        double wMin = std::numeric_limits<double>::infinity();
        double wMax = -std::numeric_limits<double>::infinity();
        for (int j = 0; j <= mesh.ny; ++j) {
            for (int i = 0; i <= mesh.nx; ++i) {
                const int n = mesh.grilleVersNoeud[j * (mesh.nx + 1) + i];
                if (n < 0) continue;
                const double w = aSol[3 * n];
                Z[j][i] = w;
                wMin = std::min(wMin, w);
                wMax = std::max(wMax, w);
            }
        }
        if (!std::isfinite(wMin) || !std::isfinite(wMax)) { wMin = 0.0; wMax = 0.0; }

        // Facteur d'amplification verticale auto-adapte : le tassement (souvent millimetrique)
        // doit rester visible a l'echelle du radier (souvent plurimetrique), sans etre fixe en
        // dur. Le facteur doit etre calibre sur l'ETENDUE du champ de tassement (wMax - wMin,
        // c'est-a-dire la profondeur reelle du bassin de deformation), et non sur sa valeur
        // absolue (max(|wMin|,|wMax|)) : un radier peut tasser presque uniformement de -1.4 mm a
        // -1.9 mm (variation reelle de seulement 0.5 mm, le "bassin"), et normaliser sur la
        // valeur absolue (1.9 mm) sous-estime alors fortement l'amplification necessaire --
        // c'est cette confusion qui aplatissait visuellement le bassin, quelle que soit par
        // ailleurs la configuration de la camera ou du ratio d'aspect de la scene 3D.
        // Ratio cible (etendue amplifiee / plus grande dimension en plan) : 0.75, soit environ
        // 1.9x le ratio precedent (0.4), pour rendre la cuvette nettement plus prononcee. Ce
        // ratio ne s'applique qu'a la copie Zexag utilisee pour dessiner la surface -- Z (les
        // valeurs reelles de tassement) n'est jamais modifie, ni utilise ailleurs que pour la
        // construction de cette geometrie d'affichage (couleur, legende, info-bulle, graduations
        // de l'axe restent toutes en valeur reelle, cf. plus bas).
        constexpr double RATIO_CIBLE_RELIEF_3D = 0.75;
        double etendueTassement = wMax - wMin;
        if (etendueTassement < 1e-12) etendueTassement = 1.0;
        const double facteurExageration =
            RATIO_CIBLE_RELIEF_3D * std::max(mesh.Lx, mesh.Ly) / etendueTassement;

        std::vector<std::vector<double>> Zexag = Z;
        for (auto& ligne : Zexag)
            for (double& v : ligne) v *= facteurExageration;

        // Graduations personnalisees de l'axe Z : la surface est dessinee avec le relief
        // amplifie (Zexag), mais chaque graduation doit malgre tout afficher la valeur REELLE
        // du tassement -- sinon l'axe Z et la barre de couleur (qui, elle, utilise Zreel)
        // affichent des echelles differentes, ce qui rend le graphique incoherent a la
        // lecture et decredibilise la comparaison avec un autre outil. On plante donc les
        // graduations aux positions amplifiees correspondantes, avec le texte en valeur reelle.
        constexpr int NB_GRADUATIONS_Z = 6;
        std::ostringstream tickValsZ, tickTextZ;
        tickValsZ << "[";
        tickTextZ << "[";
        for (int k = 0; k < NB_GRADUATIONS_Z; ++k) {
            const double t = static_cast<double>(k) / (NB_GRADUATIONS_Z - 1);
            const double valeurReelle = wMin + t * (wMax - wMin);
            if (k) { tickValsZ << ","; tickTextZ << ","; }
            tickValsZ << (valeurReelle * facteurExageration);
            tickTextZ << "'" << std::scientific << std::setprecision(3) << valeurReelle << "'";
        }
        tickValsZ << "]";
        tickTextZ << "]";

        auto ecrireVecteurJS = [](std::ostream& o, const std::vector<double>& v) {
            o << "[";
            for (std::size_t i = 0; i < v.size(); ++i) {
                if (i) o << ",";
                // NaN = position hors du domaine actif (echancrure L/U) : ecrit "null" pour
                // que Plotly laisse un trou dans la surface plutot que d'afficher un faux zero.
                if (std::isnan(v[i])) o << "null"; else o << v[i];
            }
            o << "]";
            };
        auto ecrireMatriceJS = [&](std::ostream& o, const std::vector<std::vector<double>>& M) {
            o << "[";
            for (std::size_t j = 0; j < M.size(); ++j) {
                if (j) o << ",";
                ecrireVecteurJS(o, M[j]);
            }
            o << "]";
            };

        std::ofstream f(cheminFichier);
        if (!f.is_open()) {
            std::cerr << "Erreur : impossible de creer le fichier " << cheminFichier << std::endl;
            return;
        }

        f << "<!DOCTYPE html>\n<html lang=\"fr\">\n<head>\n<meta charset=\"UTF-8\">\n";
        f << "<title>Tassement du radier (3D)</title>\n";
        f << "<script src=\"https://cdn.plot.ly/plotly-2.32.0.min.js\"></script>\n";
        f << "<style>body{font-family:Arial, sans-serif; margin:0;} #plot{width:100vw; height:100vh;}</style>\n";
        f << "</head>\n<body>\n<div id=\"plot\"></div>\n<script>\n";

        f << "const X = "; ecrireVecteurJS(f, xs); f << ";\n";
        f << "const Y = "; ecrireVecteurJS(f, ys); f << ";\n";
        f << "const Zamplifie = "; ecrireMatriceJS(f, Zexag); f << ";\n";
        f << "const Zreel = "; ecrireMatriceJS(f, Z); f << ";\n";
        f << "const data = [{\n";
        f << "  type: 'surface', x: X, y: Y, z: Zamplifie, surfacecolor: Zreel,\n";
        f << "  customdata: Zreel,\n";
        f << "  hovertemplate: 'x=%{x:.2f} m<br>y=%{y:.2f} m<br>tassement=%{customdata:.4e} m<extra></extra>',\n";
        f << "  colorscale: 'Jet', reversescale: true,\n";
        f << "  contours: {\n";
        f << "    x: { show: true, color: '#222222', width: 1, usecolormap: false },\n";
        f << "    y: { show: true, color: '#222222', width: 1, usecolormap: false }\n";
        f << "  },\n";
        f << "  colorbar: { title: { text: 'Tassement w (m)' } }\n";
        f << "}];\n";

        f << "const layout = {\n";
        f << "  title: { text: 'Tassement du radier (w) -- " << mesh.elements.size()
            << " elements actifs (grille de fond " << mesh.nx << " x " << mesh.ny << ")' },\n";
        f << "  scene: {\n";
        f << "    xaxis: { title: { text: 'X (m)' } },\n";
        f << "    yaxis: { title: { text: 'Y (m)' } },\n";
        f << "    zaxis: { title: { text: 'Tassement w (m)' },\n";
        f << "      tickvals: "; f << tickValsZ.str(); f << ",\n";
        f << "      ticktext: "; f << tickTextZ.str(); f << "\n";
        f << "    },\n";
        f << std::setprecision(6);
        // aspectmode 'data' : chaque axe est dimensionne proportionnellement a l'etendue
        // reelle des donnees affichees (y compris le relief Z deja amplifie ci-dessus), au
        // lieu d'un ratio manuel fixe qui aplatissait le bassin de tassement quelle que soit
        // la geometrie du radier.
        f << "    aspectmode: 'data',\n";
        // Camera plus rasante que la vue par defaut de Plotly (eye z abaisse) : accentue la
        // perception du creux sans aller jusqu'a une vue de profil qui masquerait la forme
        // d'ensemble du radier en plan.
        f << "    camera: { eye: { x: 1.4, y: 1.4, z: 0.65 } }\n";
        f << "  },\n";
        f << "  margin: { l: 0, r: 0, b: 0, t: 40 }\n";
        f << "};\n";

        f << "Plotly.newPlot('plot', data, layout, {responsive: true});\n";
        f << "</script>\n</body>\n</html>\n";
        f.close();

        std::cout << "\nGraphique 3D interactif du tassement exporte (HTML, a ouvrir dans un navigateur) : "
            << cheminFichier << std::endl;
    }

    // ---- Ecriture d'un vrai fichier .xlsx (archive ZIP + parties OOXML) -----
    // Contrairement a un fichier SpreadsheetML (XML) renomme en .xls, ce format
    // est le format binaire natif d'Excel : aucun avertissement de securite sur
    // la difference format/extension, aucun risque de blocage par les regles de
    // securite d'Excel/Microsoft 365. Implementation sans dependance externe :
    // les parties sont stockees dans le ZIP sans compression (methode 0), ce qui
    // evite d'avoir a implementer un algorithme de deflation.
    std::uint32_t crc32Xlsx(const std::string& data) {
        static const std::array<std::uint32_t, 256> table = [] {
            std::array<std::uint32_t, 256> t{};
            for (std::uint32_t i = 0; i < 256; ++i) {
                std::uint32_t c = i;
                for (int k = 0; k < 8; ++k) c = (c & 1u) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
                t[i] = c;
            }
            return t;
            }();
        std::uint32_t c = 0xFFFFFFFFu;
        for (unsigned char ch : data) c = table[(c ^ ch) & 0xFFu] ^ (c >> 8);
        return c ^ 0xFFFFFFFFu;
    }

    void put16Xlsx(std::string& s, unsigned v) {
        s.push_back(static_cast<char>(v & 0xFF));
        s.push_back(static_cast<char>((v >> 8) & 0xFF));
    }
    void put32Xlsx(std::string& s, std::uint32_t v) {
        s.push_back(static_cast<char>(v & 0xFF));
        s.push_back(static_cast<char>((v >> 8) & 0xFF));
        s.push_back(static_cast<char>((v >> 16) & 0xFF));
        s.push_back(static_cast<char>((v >> 24) & 0xFF));
    }

    // Tente d'ecrire l'archive ; renvoie false si le fichier ne peut pas etre
    // ouvert en ecriture (par exemple parce qu'il est deja ouvert dans Excel).
    bool ecrireXlsx(const std::string& cheminFichier,
        const std::vector<std::pair<std::string, std::string>>& parties) {
        std::ofstream out(cheminFichier, std::ios::binary | std::ios::trunc);
        if (!out.is_open()) return false;

        std::string central;
        std::uint32_t offset = 0;

        for (const auto& partie : parties) {
            const std::string& nom = partie.first;
            const std::string& contenu = partie.second;
            const std::uint32_t crc = crc32Xlsx(contenu);
            const std::uint32_t taille = static_cast<std::uint32_t>(contenu.size());

            std::string lfh;
            put32Xlsx(lfh, 0x04034b50u);
            put16Xlsx(lfh, 20);
            put16Xlsx(lfh, 0);
            put16Xlsx(lfh, 0);
            put16Xlsx(lfh, 0);
            put16Xlsx(lfh, 0x21);
            put32Xlsx(lfh, crc);
            put32Xlsx(lfh, taille);
            put32Xlsx(lfh, taille);
            put16Xlsx(lfh, static_cast<unsigned>(nom.size()));
            put16Xlsx(lfh, 0);
            lfh += nom;

            out.write(lfh.data(), static_cast<std::streamsize>(lfh.size()));
            out.write(contenu.data(), static_cast<std::streamsize>(contenu.size()));

            std::string cdh;
            put32Xlsx(cdh, 0x02014b50u);
            put16Xlsx(cdh, 20);
            put16Xlsx(cdh, 20);
            put16Xlsx(cdh, 0);
            put16Xlsx(cdh, 0);
            put16Xlsx(cdh, 0);
            put16Xlsx(cdh, 0x21);
            put32Xlsx(cdh, crc);
            put32Xlsx(cdh, taille);
            put32Xlsx(cdh, taille);
            put16Xlsx(cdh, static_cast<unsigned>(nom.size()));
            put16Xlsx(cdh, 0);
            put16Xlsx(cdh, 0);
            put16Xlsx(cdh, 0);
            put16Xlsx(cdh, 0);
            put32Xlsx(cdh, 0);
            put32Xlsx(cdh, offset);
            cdh += nom;
            central += cdh;

            offset += static_cast<std::uint32_t>(lfh.size() + contenu.size());
        }

        const std::uint32_t centralOffset = offset;
        out.write(central.data(), static_cast<std::streamsize>(central.size()));

        std::string eocd;
        put32Xlsx(eocd, 0x06054b50u);
        put16Xlsx(eocd, 0);
        put16Xlsx(eocd, 0);
        put16Xlsx(eocd, static_cast<unsigned>(parties.size()));
        put16Xlsx(eocd, static_cast<unsigned>(parties.size()));
        put32Xlsx(eocd, static_cast<std::uint32_t>(central.size()));
        put32Xlsx(eocd, centralOffset);
        put16Xlsx(eocd, 0);
        out.write(eocd.data(), static_cast<std::streamsize>(eocd.size()));

        return static_cast<bool>(out);
    }

    // Exporte tous les resultats du calcul dans un classeur .xlsx a 3 feuilles :
    // Resume, Resultats noeuds, Moments elements. Si le fichier cible est
    // verrouille (deja ouvert dans Excel), l'export bascule automatiquement sur
    // un nom de secours pour ne jamais perdre silencieusement les resultats.
    void exporterResultatsExcel(const Mesh& mesh, const PlateProps& plate,
        const std::vector<Layer>& layers, const std::vector<SurfaceLoad>& surfaceLoads,
        const std::vector<PointLoad>& pointLoads, const std::vector<LineLoad>& lineLoads,
        const SsiResult& r, const Eigen::VectorXd& Rs, const Eigen::VectorXd& ks,
        double totalCharge, const std::string& formeNom, const std::string& formeDetails,
        const std::string& cheminFichier) {

        auto xmlEchappe = [](const std::string& s) {
            std::string out;
            out.reserve(s.size());
            for (char ch : s) {
                switch (ch) {
                case '&': out += "&amp;"; break;
                case '<': out += "&lt;"; break;
                case '>': out += "&gt;"; break;
                default: out += ch;
                }
            }
            return out;
        };

        // Construit une feuille (sheetN.xml) a partir de lignes de cellules :
        // chaque cellule est du texte (t=inlineStr), un indice entier (Noeud,
        // Element : format General, pas de notation scientifique) ou une valeur
        // physique (style "s2", format scientifique a 6 chiffres significatifs).
        struct Cellule { bool estTexte; std::string texte; double valeur; bool entete; bool indice; };
        auto txt = [](std::string s, bool entete = false) {
            return Cellule{ true, std::move(s), 0.0, entete, false };
            };
        auto num = [](double v) { return Cellule{ false, "", v, false, false }; };
        auto idx = [](double v) { return Cellule{ false, "", v, false, true }; };

        // Largeur d'affichage estimee d'une cellule (pour l'ajustement auto des
        // colonnes, calcule cote C++ puisque l'ecriture est directe, sans Excel).
        auto largeurAffichee = [](const Cellule& c) -> std::size_t {
            if (c.estTexte) return c.texte.size();
            if (c.indice) return std::to_string(static_cast<long long>(c.valeur)).size();
            std::ostringstream oss;
            oss << std::scientific << std::setprecision(6) << c.valeur;
            return oss.str().size();
            };

        auto construireFeuille = [&](const std::vector<std::vector<Cellule>>& lignes) {
            // Largeurs de colonnes auto-ajustees a partir du contenu reel.
            std::size_t nbCols = 0;
            for (const auto& ligne : lignes) nbCols = std::max(nbCols, ligne.size());
            std::vector<std::size_t> largeurs(nbCols, 6);
            for (const auto& ligne : lignes) {
                for (std::size_t j = 0; j < ligne.size(); ++j)
                    largeurs[j] = std::max(largeurs[j], largeurAffichee(ligne[j]));
            }

            std::ostringstream f;
            f << std::setprecision(10);
            f << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
                "<worksheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">";
            // Quadrillage explicitement visible (comportement par defaut d'Excel,
            // rendu explicite pour ne pas dependre des reglages utilisateur).
            f << "<sheetViews><sheetView showGridLines=\"1\" workbookViewId=\"0\"/></sheetViews>";
            f << "<cols>";
            for (std::size_t j = 0; j < nbCols; ++j) {
                f << "<col min=\"" << (j + 1) << "\" max=\"" << (j + 1) << "\" width=\""
                    << (largeurs[j] + 2) << "\" customWidth=\"1\"/>";
            }
            f << "</cols>";
            f << "<sheetData>\n";
            int numeroLigne = 0;
            for (const auto& ligne : lignes) {
                ++numeroLigne;
                f << "<row r=\"" << numeroLigne << "\">";
                for (const auto& c : ligne) {
                    f << "<c";
                    if (c.entete) f << " s=\"1\"";
                    else if (!c.estTexte && !c.indice) f << " s=\"2\"";
                    if (c.estTexte) {
                        f << " t=\"inlineStr\"><is><t>" << xmlEchappe(c.texte) << "</t></is></c>";
                    }
                    else {
                        f << "><v>" << c.valeur << "</v></c>";
                    }
                }
                f << "</row>\n";
            }
            f << "</sheetData></worksheet>\n";
            return f.str();
        };

        // ---- Feuille 1 : Resume --------------------------------------------
        std::vector<std::vector<Cellule>> resume;
        resume.push_back({ txt("Resultats GeoInfraRadier") });
        resume.push_back({});
        resume.push_back({ txt("Parametre", true), txt("Valeur", true) });
        resume.push_back({ txt("Forme du radier"), txt(formeNom) });
        resume.push_back({ txt("Dimensions caracteristiques echancrure"), txt(formeDetails) });
        resume.push_back({ txt("Longueur (boite englobante) Lx (m)"), num(mesh.Lx) });
        resume.push_back({ txt("Largeur (boite englobante) Ly (m)"), num(mesh.Ly) });
        resume.push_back({ txt("Nb cellules grille de fond en X"), idx(mesh.nx) });
        resume.push_back({ txt("Nb cellules grille de fond en Y"), idx(mesh.ny) });
        resume.push_back({ txt("Nb elements actifs"), idx(static_cast<double>(mesh.elements.size())) });
        resume.push_back({ txt("Nb elements elimines (hors contour)"), idx(mesh.nbCellulesEliminees()) });
        resume.push_back({ txt("Nb noeuds actifs"), idx(mesh.nn) });
        resume.push_back({ txt("Nb ddl (ndof)"), idx(mesh.ndof) });
        resume.push_back({ txt("Module de Young beton E (kPa)"), num(plate.E) });
        resume.push_back({ txt("Coefficient de Poisson beton nu"), num(plate.nu) });
        resume.push_back({ txt("Epaisseur radier h (m)"), num(plate.h) });
        resume.push_back({});
        resume.push_back({ txt("Stratigraphie du sol", true) });
        resume.push_back({ txt("Couche", true), txt("z haut (m)", true), txt("z bas (m)", true),
                            txt("E sol (kPa)", true), txt("nu sol", true) });
        for (std::size_t i = 0; i < layers.size(); ++i) {
            const auto& L = layers[i];
            std::vector<Cellule> ligne{ idx(static_cast<double>(i + 1)), num(L.zt) };
            ligne.push_back(std::isfinite(L.zb) ? num(L.zb) : txt("infini"));
            ligne.push_back(num(L.E));
            ligne.push_back(num(L.nu));
            resume.push_back(std::move(ligne));
        }
        resume.push_back({});
        resume.push_back({ txt("Nb charges surfaciques"), idx(static_cast<double>(surfaceLoads.size())) });
        resume.push_back({ txt("Nb charges ponctuelles"), idx(static_cast<double>(pointLoads.size())) });
        resume.push_back({ txt("Nb charges lineiques"), idx(static_cast<double>(lineLoads.size())) });
        resume.push_back({});
        resume.push_back({ txt("Charge totale appliquee (kN)"), num(totalCharge) });
        resume.push_back({ txt("Somme des reactions du sol P (kN)"), num(r.P.sum()) });
        double sommeRzResume = 0.0;
        for (int i = 0; i < mesh.nn; ++i) sommeRzResume += Rs[3 * i];
        resume.push_back({ txt("Somme des efforts verticaux Rz (kN)"), num(sommeRzResume) });

        // ---- Feuille 2 : Vecteur des forces nodales F (3n x 1, 3 ddl par noeud) --
        std::vector<std::vector<Cellule>> forcesF;
        forcesF.push_back({ txt("Noeud", true), txt("Effort vert. Fz (kN)", true),
                             txt("Moment Mx (kN.m)", true), txt("Moment My (kN.m)", true) });
        for (int i = 0; i < mesh.nn; ++i) {
            forcesF.push_back({ idx(i + 1), num(r.F[3 * i]), num(r.F[3 * i + 1]),
                                 num(r.F[3 * i + 2]) });
        }

        // ---- Feuille 3 : Vecteur deplacement (3n x 1, 3 ddl par noeud) ----
        std::vector<std::vector<Cellule>> deplacements;
        deplacements.push_back({ txt("Noeud", true), txt("Fleche w (m)", true),
                                  txt("Rotation th_x (rad)", true), txt("Rotation th_y (rad)", true) });
        for (int i = 0; i < mesh.nn; ++i) {
            deplacements.push_back({ idx(i + 1), num(r.a[3 * i]), num(r.a[3 * i + 1]),
                                      num(r.a[3 * i + 2]) });
        }

        // ---- Feuille 4 : Vecteur reaction du sol Rs = Mp.p (3n x 1) -------
        std::vector<std::vector<Cellule>> reactionsRs;
        reactionsRs.push_back({ txt("Noeud", true), txt("Effort vert. Rz (kN)", true),
                                 txt("Moment Mx (kN.m)", true), txt("Moment My (kN.m)", true) });
        for (int i = 0; i < mesh.nn; ++i) {
            reactionsRs.push_back({ idx(i + 1), num(Rs[3 * i]), num(Rs[3 * i + 1]), num(Rs[3 * i + 2]) });
        }

        // ---- Feuille 5 : Vecteur pression du sol p (n x 1) -----------------
        std::vector<std::vector<Cellule>> pressionP;
        pressionP.push_back({ txt("Noeud", true), txt("Pression p (kPa)", true) });
        for (int i = 0; i < mesh.nn; ++i) pressionP.push_back({ idx(i + 1), num(r.p[i]) });

        // ---- Feuille 6 : Coefficient de reaction du sol ks = p/s (n x 1) --
        std::vector<std::vector<Cellule>> coeffKs;
        coeffKs.push_back({ txt("Noeud", true), txt("Coeff ks (kN/m3)", true) });
        for (int i = 0; i < mesh.nn; ++i) coeffKs.push_back({ idx(i + 1), num(ks[i]) });

        // ---- Feuille 7 : Moments flechissants au centre de chaque element --
        std::vector<std::vector<Cellule>> moments;
        moments.push_back({ txt("Element", true), txt("x centre (m)", true), txt("y centre (m)", true),
                             txt("Mx (kN.m par m)", true), txt("My (kN.m par m)", true),
                             txt("Mxy (kN.m par m)", true) });
        for (int e = 0; e < static_cast<int>(mesh.elements.size()); ++e) {
            double xc = 0.0, yc = 0.0;
            for (int n : mesh.elements[e]) { xc += mesh.nodes[n].x(); yc += mesh.nodes[n].y(); }
            xc /= 4.0; yc /= 4.0;
            const Eigen::Vector3d M = elementMoments(mesh, r.a, plate, e, 0.0, 0.0);
            moments.push_back({ idx(e + 1), num(xc), num(yc), num(M(0)), num(M(1)), num(M(2)) });
        }

        // ---- Assemblage des parties OOXML de l'archive .xlsx ---------------
        const std::string contentTypes =
            "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
            "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">"
            "<Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>"
            "<Default Extension=\"xml\" ContentType=\"application/xml\"/>"
            "<Override PartName=\"/xl/workbook.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml\"/>"
            "<Override PartName=\"/xl/styles.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml\"/>"
            "<Override PartName=\"/xl/worksheets/sheet1.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>"
            "<Override PartName=\"/xl/worksheets/sheet2.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>"
            "<Override PartName=\"/xl/worksheets/sheet3.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>"
            "<Override PartName=\"/xl/worksheets/sheet4.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>"
            "<Override PartName=\"/xl/worksheets/sheet5.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>"
            "<Override PartName=\"/xl/worksheets/sheet6.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>"
            "<Override PartName=\"/xl/worksheets/sheet7.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>"
            "</Types>\n";

        const std::string rootRels =
            "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
            "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
            "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" Target=\"xl/workbook.xml\"/>"
            "</Relationships>\n";

        const std::string workbookXml =
            "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
            "<workbook xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" "
            "xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\">"
            "<sheets>"
            "<sheet name=\"Resume\" sheetId=\"1\" r:id=\"rId1\"/>"
            "<sheet name=\"Forces F\" sheetId=\"2\" r:id=\"rId2\"/>"
            "<sheet name=\"Deplacements\" sheetId=\"3\" r:id=\"rId3\"/>"
            "<sheet name=\"Reactions sol Rs\" sheetId=\"4\" r:id=\"rId4\"/>"
            "<sheet name=\"Pression sol p\" sheetId=\"5\" r:id=\"rId5\"/>"
            "<sheet name=\"Coefficient ks\" sheetId=\"6\" r:id=\"rId6\"/>"
            "<sheet name=\"Moments elements\" sheetId=\"7\" r:id=\"rId7\"/>"
            "</sheets></workbook>\n";

        const std::string workbookRels =
            "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
            "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
            "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet1.xml\"/>"
            "<Relationship Id=\"rId2\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet2.xml\"/>"
            "<Relationship Id=\"rId3\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet3.xml\"/>"
            "<Relationship Id=\"rId4\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet4.xml\"/>"
            "<Relationship Id=\"rId5\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet5.xml\"/>"
            "<Relationship Id=\"rId6\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet6.xml\"/>"
            "<Relationship Id=\"rId7\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet7.xml\"/>"
            "<Relationship Id=\"rId8\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles\" Target=\"styles.xml\"/>"
            "</Relationships>\n";

        // Style 0 = normal (General). Style 1 = en-tete (gras + fond bleu clair).
        // Style 2 = valeur physique en notation scientifique a 6 decimales.
        const std::string stylesXml =
            "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
            "<styleSheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">"
            "<numFmts count=\"1\"><numFmt numFmtId=\"164\" formatCode=\"0.000000E+00\"/></numFmts>"
            "<fonts count=\"2\"><font><sz val=\"11\"/><name val=\"Calibri\"/></font>"
            "<font><b/><sz val=\"11\"/><name val=\"Calibri\"/></font></fonts>"
            "<fills count=\"2\"><fill><patternFill patternType=\"none\"/></fill>"
            "<fill><patternFill patternType=\"solid\"><fgColor rgb=\"FFD9E1F2\"/>"
            "<bgColor indexed=\"64\"/></patternFill></fill></fills>"
            "<borders count=\"1\"><border><left/><right/><top/><bottom/><diagonal/></border></borders>"
            "<cellStyleXfs count=\"1\"><xf numFmtId=\"0\" fontId=\"0\" fillId=\"0\" borderId=\"0\"/></cellStyleXfs>"
            "<cellXfs count=\"3\">"
            "<xf numFmtId=\"0\" fontId=\"0\" fillId=\"0\" borderId=\"0\" xfId=\"0\"/>"
            "<xf numFmtId=\"0\" fontId=\"1\" fillId=\"1\" borderId=\"0\" xfId=\"0\" applyFont=\"1\" applyFill=\"1\"/>"
            "<xf numFmtId=\"164\" fontId=\"0\" fillId=\"0\" borderId=\"0\" xfId=\"0\" applyNumberFormat=\"1\"/>"
            "</cellXfs></styleSheet>\n";

        const std::vector<std::pair<std::string, std::string>> parties = {
            { "[Content_Types].xml", contentTypes },
            { "_rels/.rels", rootRels },
            { "xl/workbook.xml", workbookXml },
            { "xl/_rels/workbook.xml.rels", workbookRels },
            { "xl/styles.xml", stylesXml },
            { "xl/worksheets/sheet1.xml", construireFeuille(resume) },
            { "xl/worksheets/sheet2.xml", construireFeuille(forcesF) },
            { "xl/worksheets/sheet3.xml", construireFeuille(deplacements) },
            { "xl/worksheets/sheet4.xml", construireFeuille(reactionsRs) },
            { "xl/worksheets/sheet5.xml", construireFeuille(pressionP) },
            { "xl/worksheets/sheet6.xml", construireFeuille(coeffKs) },
            { "xl/worksheets/sheet7.xml", construireFeuille(moments) },
        };

        if (ecrireXlsx(cheminFichier, parties)) {
            std::cout << "\nResultats exportes au format Excel (.xlsx, 7 feuilles : Resume, "
                "Forces F, Deplacements, Reactions sol Rs, Pression sol p, Coefficient ks, "
                "Moments elements) : " << cheminFichier << std::endl;
            return;
        }

        // Le fichier cible est probablement ouvert dans Excel (verrouille) :
        // on ecrit sous un nom de secours plutot que d'echouer silencieusement.
        for (int tentative = 2; tentative <= 9; ++tentative) {
            const std::size_t point = cheminFichier.find_last_of('.');
            const std::string secours = (point == std::string::npos)
                ? cheminFichier + "_" + std::to_string(tentative)
                : cheminFichier.substr(0, point) + "_" + std::to_string(tentative) +
                cheminFichier.substr(point);
            if (ecrireXlsx(secours, parties)) {
                std::cerr << "\nAttention : " << cheminFichier
                    << " est probablement deja ouvert dans Excel (fichier verrouille).\n"
                    << "Resultats exportes sous un nom de secours a la place : " << secours
                    << "\nFermez le classeur ouvert avant la prochaine execution pour "
                    "retrouver le nom habituel.\n" << std::endl;
                return;
            }
        }
        std::cerr << "\nErreur : impossible d'exporter les resultats vers " << cheminFichier
            << " (ni vers un nom de secours). Verifiez les droits d'ecriture du dossier.\n"
            << std::endl;
    }

    double lireDouble(const std::string& nom, double defaut) {
        while (true) {
            std::string ligne = lireLigne(nom + " [" + std::to_string(defaut) + "] : ");
            if (ligne.empty()) return defaut;
            try {
                return std::stod(ligne);
            }
            catch (...) {
                std::cout << "Valeur invalide, reessayez.\n";
            }
        }
    }

    int lireEntier(const std::string& nom, int defaut) {
        while (true) {
            std::string ligne = lireLigne(nom + " [" + std::to_string(defaut) + "] : ");
            if (ligne.empty()) return defaut;
            try {
                return std::stoi(ligne);
            }
            catch (...) {
                std::cout << "Valeur invalide, reessayez.\n";
            }
        }
    }

}  // namespace

int main() {
    // Force le vidage du tampon de sortie apres chaque insertion : sans ca,
    // les invites (qui ne se terminent pas par un saut de ligne) peuvent
    // rester bufferisees et ne jamais s'afficher avant que le programme lise
    // stdin, ce qui donne l'impression que la cellule "ne repond pas".
    std::cout.setf(std::ios::unitbuf);
    std::cout << std::fixed << std::setprecision(6);

    // ---- 1. Geometrie et maillage du radier -------------------------------
    std::cout << "===== Geometrie du radier (Entree = valeur par defaut) =====\n";
    // La forme n'a volontairement pas de valeur par defaut : l'utilisateur doit choisir
    // explicitement, une ligne vide ou une valeur hors de {1,2,3} redemande la saisie.
    std::cout << "Forme du radier : 1 = Rectangulaire, 2 = L, 3 = U\n";
    int formeChoisie = 0;
    while (true) {
        const std::string ligne = lireLigne("Choix de la forme (1, 2 ou 3) : ");
        try {
            formeChoisie = std::stoi(ligne);
            if (formeChoisie >= 1 && formeChoisie <= 3) break;
        }
        catch (...) {
            // ligne vide ou non numerique : traite comme un choix invalide ci-dessous.
        }
        std::cout << "Choix invalide : entrez 1 (Rectangulaire), 2 (L) ou 3 (U).\n";
    }

    std::string formeNom = "Rectangulaire";
    std::string formeDetails = "-";
    std::vector<Eigen::Vector2d> contour;
    double Lx, Ly;

    if (formeChoisie == 1) {
        // Radier rectangulaire classique : geometrie definie par ses deux longueurs.
        Lx = lireDouble("Longueur du radier Lx (m)", 5.0);
        Ly = lireDouble("Largeur du radier Ly (m)", 5.0);
        contour = Mesh::contourRectangle(Lx, Ly);
    }
    else {
        // Radier en L ou en U : geometrie definie par les sommets (x, y) du contour,
        // saisis dans l'ordre de parcours du polygone. Les valeurs par defaut
        // proposees (Entree seule) reproduisent une forme type de 5 m x 5 m, pour
        // pouvoir tester rapidement sans tout ressaisir.
        const bool estL = (formeChoisie == 2);
        formeNom = estL ? "L" : "U";
        const std::vector<Eigen::Vector2d> sommetsDefaut = estL
            ? Mesh::contourL(5.0, 5.0, 2.5, 2.5)
            : Mesh::contourU(5.0, 5.0, 1.25, 3.75, 2.5);

        std::cout << "-- Radier en " << formeNom << " : definissez les sommets du contour, "
            "dans l'ordre de parcours du polygone (par convention, en partant de l'origine "
            "(0,0)) --\n";
        const int nbSommetsDefaut = static_cast<int>(sommetsDefaut.size());
        int nbSommets;
        while (true) {
            nbSommets = lireEntier("Nombre de sommets du contour", nbSommetsDefaut);
            if (nbSommets >= 3) break;
            std::cout << "  Il faut au moins 3 sommets pour definir un polygone, reessayez.\n";
        }

        contour.clear();
        std::ostringstream detailsFlux;
        for (int s = 0; s < nbSommets; ++s) {
            const Eigen::Vector2d defaut = (s < nbSommetsDefaut) ? sommetsDefaut[s]
                                                                  : Eigen::Vector2d(0.0, 0.0);
            std::cout << "-- Sommet " << (s + 1) << " --\n";
            const double xs = lireDouble("  x (m)", defaut.x());
            const double ys = lireDouble("  y (m)", defaut.y());
            contour.emplace_back(xs, ys);
            if (s) detailsFlux << "; ";
            detailsFlux << "(" << xs << ", " << ys << ")";
        }
        formeDetails = detailsFlux.str();

        // La boite englobante [0,Lx] x [0,Ly] du maillage de fond est deduite des
        // sommets saisis (convention : le contour part de l'origine (0,0)).
        Lx = 0.0;
        Ly = 0.0;
        for (const auto& sommet : contour) {
            Lx = std::max(Lx, sommet.x());
            Ly = std::max(Ly, sommet.y());
        }
        if (Lx <= 0.0 || Ly <= 0.0) {
            std::cerr << "Erreur : le contour saisi ne definit pas une emprise valide "
                "(Lx et Ly doivent etre strictement positifs).\n";
            return 1;
        }
    }

    const double taille_x = lireDouble("Taille cible d'element en X (m)", 0.5);
    const double taille_y = lireDouble("Taille cible d'element en Y (m)", 0.5);
    int nx = static_cast<int>(std::lround(Lx / taille_x));
    int ny = static_cast<int>(std::lround(Ly / taille_y));
    nx = std::max(1, nx);
    ny = std::max(1, ny);

    Mesh mesh(Lx, Ly, nx, ny, contour);
    std::cout << "-> Grille de fond : " << nx << " x " << ny << " cellules ; forme '" << formeNom
        << "' -> " << mesh.elements.size() << " elements actifs (" << mesh.nbCellulesEliminees()
        << " elimines hors du contour), " << mesh.nn << " noeuds, ndof = " << mesh.ndof << "\n\n";
    // ---- 2. Materiau du radier (plaque de Kirchhoff-Love) ------------------
    std::cout << "===== Materiau du radier =====\n";
    const double E_beton = lireDouble("Module de Young du beton E (kPa)", 3.1e7);
    const double nu_beton = lireDouble("Coefficient de Poisson du beton nu", 0.25);
    const double h = lireDouble("Epaisseur du radier h (m)", 0.4);
    PlateProps plate{ E_beton, nu_beton, h };
    std::cout << "\n";

    // ---- 3. Stratigraphie du sol (demi-espace elastique multicouche) ------
    std::cout << "===== Stratigraphie du sol =====\n";
    const int nb_couches = lireEntier("Nombre de couches de sol", 1);
    std::vector<Layer> layers;
    double zTop = 0.0;
    for (int i = 0; i < nb_couches; ++i) {
        std::cout << "-- Couche " << (i + 1) << " --\n";
        const double epaisseur = lireDouble(
            "  Epaisseur compressible H (m) [0 = demi-espace infini, substratum rigide sinon]", 8.0);
        const double E_sol = lireDouble("  Module de Young E_sol (kPa)", 1.2e4);
        const double nu_sol = lireDouble("  Coefficient de Poisson nu_sol", 0.33);
        if (epaisseur <= 0.0) {
            layers.push_back({ zTop, std::numeric_limits<double>::infinity(), E_sol, nu_sol });
            std::cout << "  -> demi-espace infini sous " << zTop << " m, couches suivantes ignorees\n";
            break;
        }
        const double zBot = zTop + epaisseur;
        layers.push_back({ zTop, zBot, E_sol, nu_sol });
        zTop = zBot;
    }
    std::cout << "\n";

    // ---- 4. Chargement -------------------------------------------------
    std::cout << "===== Chargement =====\n";
    const int nb_surfaces = lireEntier("Nombre de charges surfaciques", 0);
    std::vector<SurfaceLoad> surfaceLoads;
    for (int i = 0; i < nb_surfaces; ++i) {
        std::cout << "-- Charge surfacique " << (i + 1) << " --\n";
        double xminS, xmaxS, yminS, ymaxS;
        while (true) {
            xminS = lireDouble("  Coin inferieur gauche x_min (m)", 0.0);
            yminS = lireDouble("  Coin inferieur gauche y_min (m)", 0.0);
            xmaxS = lireDouble("  Coin superieur droit x_max (m)", 0);
            ymaxS = lireDouble("  Coin superieur droit y_max (m)", 0);
            if (xmaxS > xminS && ymaxS > yminS) break;
            std::cout << "  Rectangle invalide (x_max doit etre > x_min et y_max > y_min), reessayez.\n";
        }
        const double qSurf = lireDouble("  Intensite q (kPa, negative vers le bas)", 0.0);
        surfaceLoads.push_back({ xminS, xmaxS, yminS, ymaxS, qSurf });
    }
    const int nb_charges = lireEntier("Nombre de charges ponctuelles", 1);
    std::vector<PointLoad> pointLoads;
    for (int i = 0; i < nb_charges; ++i) {
        std::cout << "-- Charge ponctuelle " << (i + 1) << " --\n";
        double xPt, yPt;
        while (true) {
            xPt = lireDouble("  Coordonnee x du point d'application (m)", 2.5);
            yPt = lireDouble("  Coordonnee y du point d'application (m)", 2.5);
            if (xPt >= 0.0 && xPt <= Lx && yPt >= 0.0 && yPt <= Ly &&
                pointDansOuSurPolygone(xPt, yPt, contour)) {
                break;
            }
            std::cout << "  Point hors du radier (hors de la boite englobante, ou dans une "
                "echancrure du L/U), reessayez.\n";
        }
        const double P = lireDouble("  Force verticale P (kN, negative vers le bas)", -200.0);
        pointLoads.push_back({ xPt, yPt, P });
    }
    const int nb_lignes = lireEntier("Nombre de charges lineiques", 1);
    std::vector<LineLoad> lineLoads;
    for (int i = 0; i < nb_lignes; ++i) {
        std::cout << "-- Charge lineique " << (i + 1) << " --\n";
        const double x1c = lireDouble("  Point de depart x1 (m)", 0.0);
        const double y1c = lireDouble("  Point de depart y1 (m)", 2.5);
        const double x2c = lireDouble("  Point d'arrivee x2 (m)", 5.0);
        const double y2c = lireDouble("  Point d'arrivee y2 (m)", 2.5);
        const double qLigne = lireDouble("  Intensite q (kN/m, negative vers le bas)", 0.0);
        lineLoads.push_back({ x1c, y1c, x2c, y2c, qLigne });
    }
    const bool afficherToutesMatrices =
        lireEntier("Afficher toutes les matrices intermediaires (1=oui, 0=non)", 1) != 0;
    std::cout << "\n";

    // ---- 5. Resolution du systeme couple sol-structure ---------------------
    SsiResult r = solveSsi(mesh, plate, layers, surfaceLoads, pointLoads, lineLoads);

    // Vecteur reaction du sol Rs = Mp * p (3n x 1 : effort vertical + moments
    // par noeud, coherent avec les 3 ddl (w, theta_x, theta_y) du vecteur a).
    const Eigen::VectorXd Rs = r.Mp * r.p;
    // Coefficient de reaction du sol ks(x,y) = p(x,y) / s(x,y) (n x 1), avec
    // securite contre la division par un tassement nul ou quasi nul.
    const Eigen::VectorXd ks = coefficientReactionSol(r.p, r.sbar);

    if (afficherToutesMatrices) {
        afficherMatrice("Matrice de rigidite de la plaque K", r.K);
        afficherMatrice("Matrice de souplesse du sol Shat", r.Shat);
        afficherMatrice("Matrice de passage Mp (ddl plaque -> forces de panneau)", r.Mp);
        afficherMatrice("Matrice de passage G (ddl plaque -> tassements de panneau)", r.G);
        std::cout << "\n";
    }

    // ---- Tableaux structures des resultats nodaux (console) ----------------
    afficherTableauVecteur3Ddl("TABLEAU 1 : Vecteur des forces nodales F (chargement applique)",
        { "Effort vert. Fz (kN)", "Moment Mx (kN.m)", "Moment My (kN.m)" }, mesh.nn, r.F);
    afficherTableauVecteur3Ddl("TABLEAU 2 : Vecteur deplacement (plaque)",
        { "Fleche w (m)", "Rotation th_x (rad)", "Rotation th_y (rad)" }, mesh.nn, r.a);
    afficherTableauVecteur3Ddl("TABLEAU 3 : Vecteur reaction du sol Rs = Mp.p",
        { "Effort vert. Rz (kN)", "Moment Mx (kN.m)", "Moment My (kN.m)" }, mesh.nn, Rs);
    afficherTableauVecteurNoeud("TABLEAU 4 : Vecteur pression du sol p", "Pression p (kPa)",
        mesh.nn, r.p);
    afficherTableauVecteurNoeud("TABLEAU 5 : Coefficient de reaction du sol ks = p/s",
        "Coeff ks (kN/m3)", mesh.nn, ks);

    double totalCharge = 0.0;
    // Approximation : suppose que chaque charge surfacique/lineique est entierement comprise
    // dans l'emprise du radier (sinon seule la portion effectivement clippee dans le maillage
    // contribue reellement a r.F, et cette somme nominale la surestimerait).
    for (const auto& sl : surfaceLoads) {
        totalCharge += sl.q * (sl.xmax - sl.xmin) * (sl.ymax - sl.ymin);
    }
    for (const auto& pl : pointLoads) totalCharge += pl.P;
    for (const auto& ll : lineLoads) {
        const double longueur = std::sqrt((ll.x2 - ll.x1) * (ll.x2 - ll.x1) +
            (ll.y2 - ll.y1) * (ll.y2 - ll.y1));
        totalCharge += ll.q * longueur;
    }
    double sommeRz = 0.0;
    for (int i = 0; i < mesh.nn; ++i) sommeRz += Rs[3 * i];
    std::cout << "\n===== Bilan d'equilibre =====\n";
    std::cout << "Somme des reactions du sol P (kN)          : " << r.P.sum() << "\n";
    std::cout << "Somme des efforts verticaux Rz (kN)        : " << sommeRz << "\n";
    std::cout << "Charge totale appliquee (kN)               : " << totalCharge << "\n";

    // ---- 6. Carte d'isovaleurs du tassement (SVG) -------------------------
    // Chemin absolu du depot : le fichier est regenere (ecrase) a chaque execution, et son
    // chemin est affiche ci-dessus/ci-dessous pour pouvoir l'ouvrir immediatement.
    const std::string cheminSvg =
        "D:/1_REDACTION_MEMOIRE_2026/radier-geoinfra/Logiciel_algorithme_07_07_26/isovaleurs_tassement_ssi.svg";
    exporterIsovaleursTassementSVG(mesh, r.a, cheminSvg);

    // ---- 7. Graphique 3D interactif du tassement (HTML) --------------------
    const std::string cheminHtml3D =
        "D:/1_REDACTION_MEMOIRE_2026/radier-geoinfra/Logiciel_algorithme_07_07_26/tassement_radier_3D.html";
    exporterGraphique3DTassementHTML(mesh, r.a, cheminHtml3D);

    // ---- 8. Export des resultats au format Excel ----------------------------
    const std::string cheminExcel =
        "D:/1_REDACTION_MEMOIRE_2026/Logiciel_Memoire/GeoInfraRadier/resultats_radier.xlsx";
    exporterResultatsExcel(mesh, plate, layers, surfaceLoads, pointLoads, lineLoads, r,
        Rs, ks, totalCharge, formeNom, formeDetails, cheminExcel);

    // Ouvre automatiquement le graphique 3D dans le navigateur par defaut : un fichier HTML
    // ne peut pas s'afficher "rendu" dans l'onglet de sortie de la cellule (contrairement au
    // SVG, qui est une image statique) -- il faut un vrai navigateur pour executer son
    // JavaScript. ShellExecute delegue a l'association de fichiers Windows (.html -> navigateur
    // par defaut), donc ceci s'ouvre a chaque execution sans action manuelle.
    HINSTANCE resultat = ShellExecuteA(nullptr, "open", cheminHtml3D.c_str(), nullptr, nullptr,
        SW_SHOWNORMAL);
    if (reinterpret_cast<INT_PTR>(resultat) <= 32) {
        std::cerr << "Attention : impossible d'ouvrir automatiquement le navigateur (code "
            << reinterpret_cast<INT_PTR>(resultat)
            << "). Ouvrez le fichier ci-dessus manuellement.\n";
    }

    return 0;
}

// ============================================================================
// API D'INTEROPERABILITE C# (DLL) -- Lot 1
// ============================================================================
// Cette section n'ajoute et ne modifie AUCUNE equation ni algorithme de calcul :
// elle expose, via une interface C plate (types simples, tableaux + compteurs,
// marshalable depuis C#/P-Invoke), le meme pipeline que celui deroule par
// main() ci-dessus (Mesh -> solveSsi -> Rs = Mp.p -> ks = p/s), en appelant
// exactement les memes fonctions (Mesh, solveSsi, coefficientReactionSol,
// elementMoments, les generateurs de contour, les exporteurs Excel/SVG/HTML).
// Reste dans ce meme fichier, compile soit dans l'executable console ci-dessus
// (fonctions alors simplement inutilisees -- aucun conflit avec main()), soit
// dans un projet DLL séparé qui inclut ce fichier source tel quel.
//
// Convention d'appel : __cdecl (defaut MSVC pour extern "C"). Cote C#, declarer
// chaque import avec [DllImport(..., CallingConvention = CallingConvention.Cdecl)].
// Codes de retour : 0 = succes, != 0 = echec (message recuperable via
// GIR_DernierMessageErreur, stocke par thread appelant).
namespace {
    thread_local std::string g_girDerniereErreur;

    // Reconstruit le maillage a partir des entrees C-ABI, en appelant exactement
    // le meme constructeur Mesh que main() (aucune logique de grille dupliquee).
    std::unique_ptr<Mesh> girConstruireMesh(int nbSommets, const double* contourXY,
        double tailleX, double tailleY) {
        std::vector<Eigen::Vector2d> contour;
        contour.reserve(nbSommets);
        for (int i = 0; i < nbSommets; ++i) contour.emplace_back(contourXY[2 * i], contourXY[2 * i + 1]);
        double Lx = 0.0, Ly = 0.0;
        for (const auto& s : contour) { Lx = std::max(Lx, s.x()); Ly = std::max(Ly, s.y()); }
        int nx = std::max(1, static_cast<int>(std::lround(Lx / tailleX)));
        int ny = std::max(1, static_cast<int>(std::lround(Ly / tailleY)));
        return std::make_unique<Mesh>(Lx, Ly, nx, ny, contour);
    }
}  // namespace

extern "C" {

    __declspec(dllexport) int GIR_DernierMessageErreur(char* tampon, int tailleTampon) {
        if (!tampon || tailleTampon <= 0) return -1;
        const int n = std::min(static_cast<int>(g_girDerniereErreur.size()), tailleTampon - 1);
        std::memcpy(tampon, g_girDerniereErreur.data(), static_cast<std::size_t>(n));
        tampon[n] = '\0';
        return 0;
    }

    // ---- Generateurs de contour : identiques a Mesh::contour*, exposes pour que
    // l'IHM C# n'ait jamais a reimplementer la geometrie des echancrures L/U.
    // xyOut doit pointer vers un buffer d'au moins 16 doubles (8 sommets max.).
    __declspec(dllexport) int GIR_ContourRectangle(double Lx, double Ly, double* xyOut,
        int* nbSommetsOut) {
        const auto c = Mesh::contourRectangle(Lx, Ly);
        for (std::size_t i = 0; i < c.size(); ++i) { xyOut[2 * i] = c[i].x(); xyOut[2 * i + 1] = c[i].y(); }
        *nbSommetsOut = static_cast<int>(c.size());
        return 0;
    }
    __declspec(dllexport) int GIR_ContourL(double Lx, double Ly, double Bx, double By,
        double* xyOut, int* nbSommetsOut) {
        if (!(Bx > 0.0 && Bx < Lx && By > 0.0 && By < Ly)) {
            g_girDerniereErreur = "Echancrure en L invalide (il faut 0 < Bx < Lx et 0 < By < Ly).";
            return 1;
        }
        const auto c = Mesh::contourL(Lx, Ly, Bx, By);
        for (std::size_t i = 0; i < c.size(); ++i) { xyOut[2 * i] = c[i].x(); xyOut[2 * i + 1] = c[i].y(); }
        *nbSommetsOut = static_cast<int>(c.size());
        return 0;
    }
    __declspec(dllexport) int GIR_ContourU(double Lx, double Ly, double x1, double x2, double By,
        double* xyOut, int* nbSommetsOut) {
        if (!(x1 >= 0.0 && x1 < x2 && x2 <= Lx && By > 0.0 && By < Ly)) {
            g_girDerniereErreur = "Echancrure en U invalide (il faut 0 <= x1 < x2 <= Lx et 0 < By < Ly).";
            return 1;
        }
        const auto c = Mesh::contourU(Lx, Ly, x1, x2, By);
        for (std::size_t i = 0; i < c.size(); ++i) { xyOut[2 * i] = c[i].x(); xyOut[2 * i + 1] = c[i].y(); }
        *nbSommetsOut = static_cast<int>(c.size());
        return 0;
    }

    // ---- Entrees / sorties du calcul complet --------------------------------
    struct GIR_Couche { double epaisseur, E, nu; };                       // m, kPa, -
    struct GIR_ChargeSurfacique { double xmin, xmax, ymin, ymax, q; };    // m, m, m, m, kPa
    struct GIR_ChargePonctuelle { double x, y, P; };                      // m, m, kN
    struct GIR_ChargeLineique { double x1, y1, x2, y2, q; };              // m, m, m, m, kN/m

    struct GIR_Entrees {
        int nbSommets;                     // >= 3, sommets du contour dans l'ordre de parcours
        const double* contourXY;           // [x0,y0,x1,y1,...] (m), longueur 2*nbSommets
        double tailleX, tailleY;           // m, taille cible d'element
        double eBeton, nuBeton, hRadier;   // kPa, -, m
        int nbCouches;
        const GIR_Couche* couches;         // epaisseur <= 0 => demi-espace infini (arrete la liste)
        int nbSurfaces;
        const GIR_ChargeSurfacique* surfaces;
        int nbPonctuelles;
        const GIR_ChargePonctuelle* ponctuelles;
        int nbLineiques;
        const GIR_ChargeLineique* lineiques;
    };

    struct GIR_Dimensions { int nn, ndof, nbElements, nx, ny; };

    // Etape 1/2 : calcule uniquement les dimensions du maillage resultant, pour
    // permettre a l'appelant (C#) d'allouer les tableaux de GIR_Resultats a la
    // bonne taille avant d'appeler GIR_Calculer.
    __declspec(dllexport) int GIR_Dimensionner(const GIR_Entrees* e, GIR_Dimensions* dimOut) {
        if (!e || !dimOut || e->nbSommets < 3) {
            g_girDerniereErreur = "Entrees invalides (contour d'au moins 3 sommets requis).";
            return 1;
        }
        try {
            const auto mesh = girConstruireMesh(e->nbSommets, e->contourXY, e->tailleX, e->tailleY);
            dimOut->nn = mesh->nn;
            dimOut->ndof = mesh->ndof;
            dimOut->nbElements = static_cast<int>(mesh->elements.size());
            dimOut->nx = mesh->nx;
            dimOut->ny = mesh->ny;
            return 0;
        }
        catch (const std::exception& ex) {
            g_girDerniereErreur = ex.what();
            return 2;
        }
    }

    // Coordonnees des noeuds et connectivite des elements du maillage, SANS resoudre le
    // systeme (pas d'appel a solveSsi) : permet a l'IHM de dessiner la grille reelle du
    // maillage (section 15/16 du cahier des charges -- affichage des elements, des noeuds et
    // des axes) des que la geometrie est saisie, avant tout calcul. Tableaux DEJA ALLOUES par
    // l'appelant (dimensionnes via GIR_Dimensionner : nn pour noeudX/noeudY, nbElements*4 pour
    // elementsNoeuds).
    struct GIR_Maillage {
        double* noeudX;
        double* noeudY;
        int* elementsNoeuds;  // 4 indices de noeud par element, dans l'ordre local du moteur
    };

    __declspec(dllexport) int GIR_ObtenirMaillage(const GIR_Entrees* e, GIR_Maillage* m) {
        if (!e || !m) { g_girDerniereErreur = "Entrees/sorties nulles."; return 1; }
        try {
            const auto mesh = girConstruireMesh(e->nbSommets, e->contourXY, e->tailleX, e->tailleY);
            for (int i = 0; i < mesh->nn; ++i) {
                m->noeudX[i] = mesh->nodes[i].x();
                m->noeudY[i] = mesh->nodes[i].y();
            }
            for (int ei = 0; ei < static_cast<int>(mesh->elements.size()); ++ei) {
                for (int k = 0; k < 4; ++k) m->elementsNoeuds[4 * ei + k] = mesh->elements[ei][k];
            }
            return 0;
        }
        catch (const std::exception& ex) {
            g_girDerniereErreur = ex.what();
            return 2;
        }
    }

    struct GIR_Resultats {
        // Tableaux (taille nn), alloues et liberes par l'appelant :
        double* noeudX; double* noeudY;
        double* w; double* thetaX; double* thetaY;             // deplacement (a)
        double* fz; double* fMx; double* fMy;                  // forces nodales (F)
        double* rz; double* rMx; double* rMy;                  // reaction du sol (Rs = Mp.p)
        double* pression;                                       // p
        double* ks;                                             // coefficient de reaction (p/s)
        // Tableaux (taille nbElements) :
        double* elemXc; double* elemYc;
        double* elemMx; double* elemMy; double* elemMxy;
        // Scalaires (bilan d'equilibre) :
        double sommeReactionP;
        double sommeRz;
    };

    // Etape 2/2 : execute le calcul complet et remplit les tableaux DEJA ALLOUES
    // par l'appelant (dimensionnes via GIR_Dimensionner). Appelle exactement les
    // memes fonctions que main() : aucune equation n'est reimplementee ici.
    __declspec(dllexport) int GIR_Calculer(const GIR_Entrees* e, GIR_Resultats* r) {
        if (!e || !r) { g_girDerniereErreur = "Entrees/sorties nulles."; return 1; }
        try {
            const auto mesh = girConstruireMesh(e->nbSommets, e->contourXY, e->tailleX, e->tailleY);
            const PlateProps plate{ e->eBeton, e->nuBeton, e->hRadier };

            std::vector<Layer> layers;
            double zTop = 0.0;
            for (int i = 0; i < e->nbCouches; ++i) {
                const auto& c = e->couches[i];
                if (c.epaisseur <= 0.0) {
                    layers.push_back({ zTop, std::numeric_limits<double>::infinity(), c.E, c.nu });
                    break;
                }
                const double zBot = zTop + c.epaisseur;
                layers.push_back({ zTop, zBot, c.E, c.nu });
                zTop = zBot;
            }

            std::vector<SurfaceLoad> surfaceLoads;
            for (int i = 0; i < e->nbSurfaces; ++i) {
                const auto& s = e->surfaces[i];
                surfaceLoads.push_back({ s.xmin, s.xmax, s.ymin, s.ymax, s.q });
            }
            std::vector<PointLoad> pointLoads;
            for (int i = 0; i < e->nbPonctuelles; ++i) {
                const auto& p = e->ponctuelles[i];
                pointLoads.push_back({ p.x, p.y, p.P });
            }
            std::vector<LineLoad> lineLoads;
            for (int i = 0; i < e->nbLineiques; ++i) {
                const auto& l = e->lineiques[i];
                lineLoads.push_back({ l.x1, l.y1, l.x2, l.y2, l.q });
            }

            const SsiResult res = solveSsi(*mesh, plate, layers, surfaceLoads, pointLoads, lineLoads);
            const Eigen::VectorXd Rs = res.Mp * res.p;
            const Eigen::VectorXd ks = coefficientReactionSol(res.p, res.sbar);

            for (int i = 0; i < mesh->nn; ++i) {
                r->noeudX[i] = mesh->nodes[i].x();
                r->noeudY[i] = mesh->nodes[i].y();
                r->w[i] = res.a[3 * i];
                r->thetaX[i] = res.a[3 * i + 1];
                r->thetaY[i] = res.a[3 * i + 2];
                r->fz[i] = res.F[3 * i];
                r->fMx[i] = res.F[3 * i + 1];
                r->fMy[i] = res.F[3 * i + 2];
                r->rz[i] = Rs[3 * i];
                r->rMx[i] = Rs[3 * i + 1];
                r->rMy[i] = Rs[3 * i + 2];
                r->pression[i] = res.p[i];
                r->ks[i] = ks[i];
            }
            for (int ei = 0; ei < static_cast<int>(mesh->elements.size()); ++ei) {
                double xc = 0.0, yc = 0.0;
                for (int n : mesh->elements[ei]) { xc += mesh->nodes[n].x(); yc += mesh->nodes[n].y(); }
                xc /= 4.0; yc /= 4.0;
                const Eigen::Vector3d M = elementMoments(*mesh, res.a, plate, ei, 0.0, 0.0);
                r->elemXc[ei] = xc; r->elemYc[ei] = yc;
                r->elemMx[ei] = M(0); r->elemMy[ei] = M(1); r->elemMxy[ei] = M(2);
            }
            r->sommeReactionP = res.P.sum();
            double sommeRz = 0.0;
            for (int i = 0; i < mesh->nn; ++i) sommeRz += Rs[3 * i];
            r->sommeRz = sommeRz;
            return 0;
        }
        catch (const std::exception& ex) {
            g_girDerniereErreur = ex.what();
            return 2;
        }
    }

    // Rejoue le meme calcul et exporte, en plus, les 3 fichiers habituels
    // (xlsx/svg/html) aux chemins fournis par l'appelant, au lieu des chemins
    // fixes utilises par le mode console. Passer nullptr pour omettre un export.
    __declspec(dllexport) int GIR_Exporter(const GIR_Entrees* e, const char* cheminXlsx,
        const char* cheminSvg, const char* cheminHtml) {
        if (!e) { g_girDerniereErreur = "Entrees nulles."; return 1; }
        try {
            const auto mesh = girConstruireMesh(e->nbSommets, e->contourXY, e->tailleX, e->tailleY);
            const PlateProps plate{ e->eBeton, e->nuBeton, e->hRadier };

            std::vector<Layer> layers;
            double zTop = 0.0;
            for (int i = 0; i < e->nbCouches; ++i) {
                const auto& c = e->couches[i];
                if (c.epaisseur <= 0.0) {
                    layers.push_back({ zTop, std::numeric_limits<double>::infinity(), c.E, c.nu });
                    break;
                }
                const double zBot = zTop + c.epaisseur;
                layers.push_back({ zTop, zBot, c.E, c.nu });
                zTop = zBot;
            }
            std::vector<SurfaceLoad> surfaceLoads;
            for (int i = 0; i < e->nbSurfaces; ++i) {
                const auto& s = e->surfaces[i];
                surfaceLoads.push_back({ s.xmin, s.xmax, s.ymin, s.ymax, s.q });
            }
            std::vector<PointLoad> pointLoads;
            for (int i = 0; i < e->nbPonctuelles; ++i) {
                const auto& p = e->ponctuelles[i];
                pointLoads.push_back({ p.x, p.y, p.P });
            }
            std::vector<LineLoad> lineLoads;
            for (int i = 0; i < e->nbLineiques; ++i) {
                const auto& l = e->lineiques[i];
                lineLoads.push_back({ l.x1, l.y1, l.x2, l.y2, l.q });
            }

            const SsiResult res = solveSsi(*mesh, plate, layers, surfaceLoads, pointLoads, lineLoads);
            const Eigen::VectorXd Rs = res.Mp * res.p;
            const Eigen::VectorXd ks = coefficientReactionSol(res.p, res.sbar);

            double totalCharge = 0.0;
            for (const auto& s : surfaceLoads) totalCharge += s.q * (s.xmax - s.xmin) * (s.ymax - s.ymin);
            for (const auto& p : pointLoads) totalCharge += p.P;
            for (const auto& l : lineLoads) {
                const double longueur = std::sqrt((l.x2 - l.x1) * (l.x2 - l.x1) + (l.y2 - l.y1) * (l.y2 - l.y1));
                totalCharge += l.q * longueur;
            }

            if (cheminSvg) exporterIsovaleursTassementSVG(*mesh, res.a, cheminSvg);
            if (cheminHtml) exporterGraphique3DTassementHTML(*mesh, res.a, cheminHtml);
            if (cheminXlsx) {
                exporterResultatsExcel(*mesh, plate, layers, surfaceLoads, pointLoads, lineLoads, res,
                    Rs, ks, totalCharge, "Personnalisee (IHM)", "-", cheminXlsx);
            }
            return 0;
        }
        catch (const std::exception& ex) {
            g_girDerniereErreur = ex.what();
            return 2;
        }
    }

}  // extern "C"
