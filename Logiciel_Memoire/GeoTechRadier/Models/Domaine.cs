namespace GeoTechRadier.Models;

/// <summary>Forme du contour du radier, telle que geree par Mesh::contourRectangle/L/U dans le moteur.</summary>
public enum FormeRadier
{
    Rectangulaire,
    L,
    U
}

/// <summary>Couche de sol (m, kPa, -). Epaisseur &lt;= 0 signifie demi-espace infini (arrete la stratigraphie).</summary>
public class CoucheSol
{
    public double Epaisseur { get; set; } = 8.0;
    public double E { get; set; } = 12000.0;
    public double Nu { get; set; } = 0.33;
    public bool DemiEspaceInfini { get; set; }
}

public class ChargeSurfacique
{
    public double Xmin { get; set; }
    public double Xmax { get; set; }
    public double Ymin { get; set; }
    public double Ymax { get; set; }
    public double Q { get; set; }
}

public class ChargePonctuelle
{
    public double X { get; set; }
    public double Y { get; set; }
    public double P { get; set; } = -200.0;
}

public class ChargeLineique
{
    public double X1 { get; set; }
    public double Y1 { get; set; }
    public double X2 { get; set; }
    public double Y2 { get; set; }
    public double Q { get; set; }
}

public class Sommet
{
    public double X { get; set; }
    public double Y { get; set; }
}

/// <summary>Resultat pour un noeud du maillage -- toutes les valeurs proviennent du moteur C++
/// (GIR_Calculer), aucune n'est calculee cote C#.</summary>
public class NoeudResultat
{
    public int Numero { get; set; }
    public double X { get; set; }
    public double Y { get; set; }
    public double W { get; set; }
    public double ThetaX { get; set; }
    public double ThetaY { get; set; }
    public double Fz { get; set; }
    public double FMx { get; set; }
    public double FMy { get; set; }
    public double Rz { get; set; }
    public double RMx { get; set; }
    public double RMy { get; set; }
    public double Pression { get; set; }
    public double Ks { get; set; }
}

/// <summary>Resultat pour un element du maillage (moments flechissants au centre).</summary>
public class ElementResultat
{
    public int Numero { get; set; }
    public double Xc { get; set; }
    public double Yc { get; set; }
    public double Mx { get; set; }
    public double My { get; set; }
    public double Mxy { get; set; }
}

/// <summary>Geometrie complete du maillage (noeuds + connectivite des elements), obtenue via
/// GIR_ObtenirMaillage -- sans resolution du systeme. Utilisee pour dessiner reellement la
/// grille du maillage dans le module Radier &amp; Maillage (section 15/16 du cahier des
/// charges), et non plus seulement le contour exterieur.</summary>
public class MaillageGeometrie
{
    public List<Sommet> Noeuds { get; } = new();
    /// <summary>Chaque element sous forme de 4 sommets (coin) dans l'ordre local du moteur.</summary>
    public List<Sommet[]> Elements { get; } = new();
}

/// <summary>Dimensions du maillage genere par le moteur (section 8/16 du cahier des charges),
/// obtenues via GIR_Dimensionner seul -- disponible immediatement apres la saisie de la
/// geometrie, avant meme de lancer un calcul complet.</summary>
public class MaillageInfo
{
    public int Nn { get; set; }
    public int Ndof { get; set; }
    public int NbElements { get; set; }
    public int Nx { get; set; }
    public int Ny { get; set; }
}

/// <summary>Bilan et resultats complets d'un calcul, tels que renvoyes par GIR_Calculer.</summary>
public class ResultatsCalcul
{
    public int Nn { get; set; }
    public int Ndof { get; set; }
    public int NbElements { get; set; }
    public int Nx { get; set; }
    public int Ny { get; set; }
    public double SommeReactionP { get; set; }
    public double SommeRz { get; set; }
    public List<NoeudResultat> Noeuds { get; } = new();
    public List<ElementResultat> Elements { get; } = new();
}

/// <summary>Synthese des valeurs extremes (section 23/32 du cahier des charges). Ce n'est PAS
/// un resultat calcule par le moteur : c'est un simple post-traitement (min/max/distance) sur
/// les tableaux deja renvoyes par GIR_Calculer -- aucune nouvelle grandeur mecanique n'est
/// introduite, uniquement une lecture agregee des valeurs existantes.</summary>
public class ValeursExtremes
{
    public NoeudResultat? TassementMax { get; set; }
    public NoeudResultat? TassementMin { get; set; }
    public NoeudResultat? PressionMax { get; set; }
    public NoeudResultat? PressionMin { get; set; }
    public NoeudResultat? ReactionMax { get; set; }
    public double ReactionTotale { get; set; }

    /// <summary>Distorsion angulaire omega = |Smax - Smin| / d (section 21 du cahier des
    /// charges), avec d = distance entre les noeuds de tassement extreme.</summary>
    public double DistorsionMax { get; set; }
    public double DistanceCorrespondante { get; set; }

    public static ValeursExtremes Calculer(ResultatsCalcul r)
    {
        var v = new ValeursExtremes();
        if (r.Noeuds.Count == 0) return v;

        v.TassementMax = r.Noeuds.MaxBy(n => n.W);
        v.TassementMin = r.Noeuds.MinBy(n => n.W);
        v.PressionMax = r.Noeuds.MaxBy(n => n.Pression);
        v.PressionMin = r.Noeuds.MinBy(n => n.Pression);
        v.ReactionMax = r.Noeuds.MaxBy(n => Math.Abs(n.Rz));
        v.ReactionTotale = r.Noeuds.Sum(n => n.Rz);

        if (v.TassementMax != null && v.TassementMin != null)
        {
            var dx = v.TassementMax.X - v.TassementMin.X;
            var dy = v.TassementMax.Y - v.TassementMin.Y;
            v.DistanceCorrespondante = Math.Sqrt(dx * dx + dy * dy);
            var deltaS = Math.Abs(v.TassementMax.W - v.TassementMin.W);
            v.DistorsionMax = v.DistanceCorrespondante > 1e-9 ? deltaS / v.DistanceCorrespondante : 0.0;
        }
        return v;
    }
}
