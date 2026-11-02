using System.Runtime.InteropServices;
using System.Text;
using GeoTechRadier.Interop;
using GeoTechRadier.Models;

namespace GeoTechRadier.Services;

/// <summary>Exception levee lorsque le moteur natif (GeoInfraRadierCalcul.dll) renvoie un code
/// d'erreur : le message provient de GIR_DernierMessageErreur, donc directement du moteur C++.</summary>
public class MoteurCalculException : Exception
{
    public MoteurCalculException(string message) : base(message) { }
}

/// <summary>Donnees d'entree completes d'un calcul, cote C# (avant conversion vers la structure
/// C-ABI GIR_Entrees). Correspond exactement a la sequence de saisie de main() dans
/// GeoInfraRadier.cpp : forme -> contour -> maillage -> materiau -> sol -> chargements.</summary>
public class DonneesCalcul
{
    public FormeRadier Forme { get; set; } = FormeRadier.Rectangulaire;

    // Rectangulaire :
    public double Lx { get; set; } = 5.0;
    public double Ly { get; set; } = 5.0;
    // L :
    public double Bx { get; set; } = 2.5;
    public double By { get; set; } = 2.5;
    // U :
    public double X1 { get; set; } = 1.25;
    public double X2 { get; set; } = 3.75;
    public double ByU { get; set; } = 2.5;

    public double TailleX { get; set; } = 0.5;
    public double TailleY { get; set; } = 0.5;

    public double EBeton { get; set; } = 3.1e7;
    public double NuBeton { get; set; } = 0.25;
    public double HRadier { get; set; } = 0.4;

    public List<CoucheSol> Couches { get; } = new() { new CoucheSol() };
    public List<ChargeSurfacique> Surfaces { get; } = new();
    public List<ChargePonctuelle> Ponctuelles { get; } = new() { new ChargePonctuelle { X = 2.5, Y = 2.5 } };
    public List<ChargeLineique> Lineiques { get; } = new();
}

/// <summary>Point d'entree unique vers le moteur natif. Ne contient AUCUN calcul mecanique :
/// uniquement la conversion des donnees C# vers le contrat C-ABI, le pinning memoire necessaire
/// au marshalling P/Invoke, et la traduction des codes d'erreur natifs.</summary>
public static class RadierCalculService
{
    /// <summary>Construit le contour (liste de sommets) exactement comme le ferait le moteur,
    /// en appelant les generateurs natifs Mesh::contourRectangle/L/U -- jamais reimplementes
    /// cote C#, pour garantir zero divergence geometrique avec GeoInfraRadier.cpp.</summary>
    public static List<Sommet> ConstruireContour(DonneesCalcul d)
    {
        var buffer = new double[16];
        int rc; int nbSommets;
        switch (d.Forme)
        {
            case FormeRadier.Rectangulaire:
                rc = NativeMethods.GIR_ContourRectangle(d.Lx, d.Ly, buffer, out nbSommets);
                break;
            case FormeRadier.L:
                rc = NativeMethods.GIR_ContourL(d.Lx, d.Ly, d.Bx, d.By, buffer, out nbSommets);
                break;
            default:
                rc = NativeMethods.GIR_ContourU(d.Lx, d.Ly, d.X1, d.X2, d.ByU, buffer, out nbSommets);
                break;
        }
        if (rc != 0) throw new MoteurCalculException(DernierMessageErreur());

        var sommets = new List<Sommet>(nbSommets);
        for (int i = 0; i < nbSommets; i++)
            sommets.Add(new Sommet { X = buffer[2 * i], Y = buffer[2 * i + 1] });
        return sommets;
    }

    /// <summary>Calcule uniquement les dimensions du maillage (GIR_Dimensionner), sans resoudre
    /// le systeme -- utilise pour rafraichir les indicateurs du module Radier &amp; Maillage
    /// (section 16) des que la geometrie change, sans attendre un calcul complet.</summary>
    public static MaillageInfo Dimensionner(DonneesCalcul d)
    {
        var sommets = ConstruireContour(d);
        var contourXY = new double[sommets.Count * 2];
        for (int i = 0; i < sommets.Count; i++)
        {
            contourXY[2 * i] = sommets[i].X;
            contourXY[2 * i + 1] = sommets[i].Y;
        }
        var handle = GCHandle.Alloc(contourXY, GCHandleType.Pinned);
        try
        {
            var entrees = new GIR_Entrees
            {
                NbSommets = sommets.Count,
                ContourXY = handle.AddrOfPinnedObject(),
                TailleX = d.TailleX,
                TailleY = d.TailleY,
            };
            if (NativeMethods.GIR_Dimensionner(ref entrees, out var dim) != 0)
                throw new MoteurCalculException(DernierMessageErreur());
            return new MaillageInfo { Nn = dim.Nn, Ndof = dim.Ndof, NbElements = dim.NbElements, Nx = dim.Nx, Ny = dim.Ny };
        }
        finally
        {
            handle.Free();
        }
    }

    /// <summary>Recupere la geometrie complete du maillage (noeuds + connectivite des elements)
    /// sans resoudre le systeme, pour un affichage graphique reel de la grille (section 15/16),
    /// et non plus seulement du contour exterieur.</summary>
    public static MaillageGeometrie ObtenirMaillageGeometrie(DonneesCalcul d)
    {
        var sommets = ConstruireContour(d);
        var contourXY = new double[sommets.Count * 2];
        for (int i = 0; i < sommets.Count; i++)
        {
            contourXY[2 * i] = sommets[i].X;
            contourXY[2 * i + 1] = sommets[i].Y;
        }

        var handles = new List<GCHandle>();
        try
        {
            IntPtr Pin<T>(T[] arr) where T : struct
            {
                if (arr.Length == 0) return IntPtr.Zero;
                var h = GCHandle.Alloc(arr, GCHandleType.Pinned);
                handles.Add(h);
                return h.AddrOfPinnedObject();
            }

            var entrees = new GIR_Entrees
            {
                NbSommets = sommets.Count,
                ContourXY = Pin(contourXY),
                TailleX = d.TailleX,
                TailleY = d.TailleY,
            };

            if (NativeMethods.GIR_Dimensionner(ref entrees, out var dim) != 0)
                throw new MoteurCalculException(DernierMessageErreur());

            var noeudX = new double[dim.Nn];
            var noeudY = new double[dim.Nn];
            var elementsNoeuds = new int[dim.NbElements * 4];

            var maillage = new GIR_Maillage
            {
                NoeudX = Pin(noeudX),
                NoeudY = Pin(noeudY),
                ElementsNoeuds = Pin(elementsNoeuds),
            };

            if (NativeMethods.GIR_ObtenirMaillage(ref entrees, ref maillage) != 0)
                throw new MoteurCalculException(DernierMessageErreur());

            var geometrie = new MaillageGeometrie();
            for (int i = 0; i < dim.Nn; i++)
                geometrie.Noeuds.Add(new Sommet { X = noeudX[i], Y = noeudY[i] });
            for (int ei = 0; ei < dim.NbElements; ei++)
            {
                var quad = new Sommet[4];
                for (int k = 0; k < 4; k++)
                {
                    int idx = elementsNoeuds[4 * ei + k];
                    quad[k] = new Sommet { X = noeudX[idx], Y = noeudY[idx] };
                }
                geometrie.Elements.Add(quad);
            }
            return geometrie;
        }
        finally
        {
            foreach (var h in handles) h.Free();
        }
    }

    /// <summary>Lance le calcul complet (dimensionnement + resolution) et renvoie tous les
    /// resultats nodaux/elementaires. Chaque valeur provient d'un tableau rempli par
    /// GIR_Calculer -- aucune n'est calculee ou approximee cote C#.</summary>
    /// <summary>Etapes reelles franchies par un appel a Calculer, dans l'ordre ou elles
    /// surviennent effectivement (deux appels natifs distincts) -- utilise par le module
    /// "Calcul" de l'IHM pour un suivi honnete (pas de sous-etapes fictives que le contrat
    /// GIR_* actuel ne permet pas d'observer depuis l'exterieur du moteur).</summary>
    public enum EtapeNative { PreparationMaillage, ResolutionSysteme }

    public static ResultatsCalcul Calculer(DonneesCalcul d, Action<EtapeNative>? surEtape = null)
    {
        var sommets = ConstruireContour(d);
        var contourXY = new double[sommets.Count * 2];
        for (int i = 0; i < sommets.Count; i++)
        {
            contourXY[2 * i] = sommets[i].X;
            contourXY[2 * i + 1] = sommets[i].Y;
        }

        var couches = d.Couches.Select(c => new GIR_Couche
        {
            Epaisseur = c.DemiEspaceInfini ? 0.0 : c.Epaisseur,
            E = c.E,
            Nu = c.Nu
        }).ToArray();
        var surfaces = d.Surfaces.Select(s => new GIR_ChargeSurfacique
        { Xmin = s.Xmin, Xmax = s.Xmax, Ymin = s.Ymin, Ymax = s.Ymax, Q = s.Q }).ToArray();
        var ponctuelles = d.Ponctuelles.Select(p => new GIR_ChargePonctuelle
        { X = p.X, Y = p.Y, P = p.P }).ToArray();
        var lineiques = d.Lineiques.Select(l => new GIR_ChargeLineique
        { X1 = l.X1, Y1 = l.Y1, X2 = l.X2, Y2 = l.Y2, Q = l.Q }).ToArray();

        // GCHandle.Pinned : empeche le GC de deplacer les tableaux managés pendant que le
        // moteur natif detient un pointeur brut vers leur contenu (duree de vie du bloc try).
        var handles = new List<GCHandle>();
        try
        {
            IntPtr Pin<T>(T[] arr) where T : struct
            {
                if (arr.Length == 0) return IntPtr.Zero;
                var h = GCHandle.Alloc(arr, GCHandleType.Pinned);
                handles.Add(h);
                return h.AddrOfPinnedObject();
            }

            var entrees = new GIR_Entrees
            {
                NbSommets = sommets.Count,
                ContourXY = Pin(contourXY),
                TailleX = d.TailleX,
                TailleY = d.TailleY,
                EBeton = d.EBeton,
                NuBeton = d.NuBeton,
                HRadier = d.HRadier,
                NbCouches = couches.Length,
                Couches = Pin(couches),
                NbSurfaces = surfaces.Length,
                Surfaces = Pin(surfaces),
                NbPonctuelles = ponctuelles.Length,
                Ponctuelles = Pin(ponctuelles),
                NbLineiques = lineiques.Length,
                Lineiques = Pin(lineiques),
            };

            surEtape?.Invoke(EtapeNative.PreparationMaillage);
            if (NativeMethods.GIR_Dimensionner(ref entrees, out var dim) != 0)
                throw new MoteurCalculException(DernierMessageErreur());

            var noeudX = new double[dim.Nn]; var noeudY = new double[dim.Nn];
            var w = new double[dim.Nn]; var thetaX = new double[dim.Nn]; var thetaY = new double[dim.Nn];
            var fz = new double[dim.Nn]; var fMx = new double[dim.Nn]; var fMy = new double[dim.Nn];
            var rz = new double[dim.Nn]; var rMx = new double[dim.Nn]; var rMy = new double[dim.Nn];
            var pression = new double[dim.Nn]; var ks = new double[dim.Nn];
            var elemXc = new double[dim.NbElements]; var elemYc = new double[dim.NbElements];
            var elemMx = new double[dim.NbElements]; var elemMy = new double[dim.NbElements];
            var elemMxy = new double[dim.NbElements];

            var resultats = new GIR_Resultats
            {
                NoeudX = Pin(noeudX),
                NoeudY = Pin(noeudY),
                W = Pin(w),
                ThetaX = Pin(thetaX),
                ThetaY = Pin(thetaY),
                Fz = Pin(fz),
                FMx = Pin(fMx),
                FMy = Pin(fMy),
                Rz = Pin(rz),
                RMx = Pin(rMx),
                RMy = Pin(rMy),
                Pression = Pin(pression),
                Ks = Pin(ks),
                ElemXc = Pin(elemXc),
                ElemYc = Pin(elemYc),
                ElemMx = Pin(elemMx),
                ElemMy = Pin(elemMy),
                ElemMxy = Pin(elemMxy),
            };

            surEtape?.Invoke(EtapeNative.ResolutionSysteme);
            if (NativeMethods.GIR_Calculer(ref entrees, ref resultats) != 0)
                throw new MoteurCalculException(DernierMessageErreur());

            var r = new ResultatsCalcul
            {
                Nn = dim.Nn,
                Ndof = dim.Ndof,
                NbElements = dim.NbElements,
                Nx = dim.Nx,
                Ny = dim.Ny,
                SommeReactionP = resultats.SommeReactionP,
                SommeRz = resultats.SommeRz,
            };
            for (int i = 0; i < dim.Nn; i++)
            {
                r.Noeuds.Add(new NoeudResultat
                {
                    Numero = i + 1,
                    X = noeudX[i],
                    Y = noeudY[i],
                    W = w[i],
                    ThetaX = thetaX[i],
                    ThetaY = thetaY[i],
                    Fz = fz[i],
                    FMx = fMx[i],
                    FMy = fMy[i],
                    Rz = rz[i],
                    RMx = rMx[i],
                    RMy = rMy[i],
                    Pression = pression[i],
                    Ks = ks[i],
                });
            }
            for (int i = 0; i < dim.NbElements; i++)
            {
                r.Elements.Add(new ElementResultat
                {
                    Numero = i + 1,
                    Xc = elemXc[i],
                    Yc = elemYc[i],
                    Mx = elemMx[i],
                    My = elemMy[i],
                    Mxy = elemMxy[i],
                });
            }
            return r;
        }
        finally
        {
            foreach (var h in handles) h.Free();
        }
    }

    /// <summary>Rejoue le calcul et exporte, en plus, les fichiers habituels (xlsx/svg/html) du
    /// moteur, aux chemins fournis (au lieu des chemins fixes du mode console).</summary>
    public static void Exporter(DonneesCalcul d, string? cheminXlsx, string? cheminSvg, string? cheminHtml)
    {
        var sommets = ConstruireContour(d);
        var contourXY = new double[sommets.Count * 2];
        for (int i = 0; i < sommets.Count; i++)
        {
            contourXY[2 * i] = sommets[i].X;
            contourXY[2 * i + 1] = sommets[i].Y;
        }
        var couches = d.Couches.Select(c => new GIR_Couche
        { Epaisseur = c.DemiEspaceInfini ? 0.0 : c.Epaisseur, E = c.E, Nu = c.Nu }).ToArray();
        var surfaces = d.Surfaces.Select(s => new GIR_ChargeSurfacique
        { Xmin = s.Xmin, Xmax = s.Xmax, Ymin = s.Ymin, Ymax = s.Ymax, Q = s.Q }).ToArray();
        var ponctuelles = d.Ponctuelles.Select(p => new GIR_ChargePonctuelle
        { X = p.X, Y = p.Y, P = p.P }).ToArray();
        var lineiques = d.Lineiques.Select(l => new GIR_ChargeLineique
        { X1 = l.X1, Y1 = l.Y1, X2 = l.X2, Y2 = l.Y2, Q = l.Q }).ToArray();

        var handles = new List<GCHandle>();
        try
        {
            IntPtr Pin<T>(T[] arr) where T : struct
            {
                if (arr.Length == 0) return IntPtr.Zero;
                var h = GCHandle.Alloc(arr, GCHandleType.Pinned);
                handles.Add(h);
                return h.AddrOfPinnedObject();
            }

            var entrees = new GIR_Entrees
            {
                NbSommets = sommets.Count,
                ContourXY = Pin(contourXY),
                TailleX = d.TailleX,
                TailleY = d.TailleY,
                EBeton = d.EBeton,
                NuBeton = d.NuBeton,
                HRadier = d.HRadier,
                NbCouches = couches.Length,
                Couches = Pin(couches),
                NbSurfaces = surfaces.Length,
                Surfaces = Pin(surfaces),
                NbPonctuelles = ponctuelles.Length,
                Ponctuelles = Pin(ponctuelles),
                NbLineiques = lineiques.Length,
                Lineiques = Pin(lineiques),
            };

            if (NativeMethods.GIR_Exporter(ref entrees, cheminXlsx, cheminSvg, cheminHtml) != 0)
                throw new MoteurCalculException(DernierMessageErreur());
        }
        finally
        {
            foreach (var h in handles) h.Free();
        }
    }

    private static string DernierMessageErreur()
    {
        var tampon = new byte[1024];
        NativeMethods.GIR_DernierMessageErreur(tampon, tampon.Length);
        int fin = Array.IndexOf(tampon, (byte)0);
        if (fin < 0) fin = tampon.Length;
        return Encoding.ASCII.GetString(tampon, 0, fin);
    }
}
