using System.Runtime.InteropServices;

namespace GeoTechRadier.Interop;

// Ces structures reproduisent EXACTEMENT (memes champs, meme ordre) les structures C-ABI
// exportees depuis GeoInfraRadier.cpp (bloc "API D'INTEROPERABILITE C# (DLL)", fin de fichier).
// Aucune equation ni logique de maillage n'est reimplementee ici : ce fichier ne fait que
// decrire la forme binaire du contrat pour permettre au CLR de marshaler les appels.

[StructLayout(LayoutKind.Sequential)]
public struct GIR_Couche
{
    public double Epaisseur, E, Nu;
}

[StructLayout(LayoutKind.Sequential)]
public struct GIR_ChargeSurfacique
{
    public double Xmin, Xmax, Ymin, Ymax, Q;
}

[StructLayout(LayoutKind.Sequential)]
public struct GIR_ChargePonctuelle
{
    public double X, Y, P;
}

[StructLayout(LayoutKind.Sequential)]
public struct GIR_ChargeLineique
{
    public double X1, Y1, X2, Y2, Q;
}

[StructLayout(LayoutKind.Sequential)]
public struct GIR_Entrees
{
    public int NbSommets;
    public IntPtr ContourXY;          // double* -- [x0,y0,x1,y1,...]
    public double TailleX, TailleY;
    public double EBeton, NuBeton, HRadier;
    public int NbCouches;
    public IntPtr Couches;             // GIR_Couche*
    public int NbSurfaces;
    public IntPtr Surfaces;            // GIR_ChargeSurfacique*
    public int NbPonctuelles;
    public IntPtr Ponctuelles;         // GIR_ChargePonctuelle*
    public int NbLineiques;
    public IntPtr Lineiques;           // GIR_ChargeLineique*
}

[StructLayout(LayoutKind.Sequential)]
public struct GIR_Dimensions
{
    public int Nn, Ndof, NbElements, Nx, Ny;
}

[StructLayout(LayoutKind.Sequential)]
public struct GIR_Maillage
{
    public IntPtr NoeudX, NoeudY;      // double* , taille Nn
    public IntPtr ElementsNoeuds;      // int* , taille NbElements*4
}

[StructLayout(LayoutKind.Sequential)]
public struct GIR_Resultats
{
    public IntPtr NoeudX, NoeudY;
    public IntPtr W, ThetaX, ThetaY;
    public IntPtr Fz, FMx, FMy;
    public IntPtr Rz, RMx, RMy;
    public IntPtr Pression, Ks;
    public IntPtr ElemXc, ElemYc, ElemMx, ElemMy, ElemMxy;
    public double SommeReactionP, SommeRz;
}

internal static class NativeMethods
{
    // Nom de la DLL native : GeoInfraRadierCalcul.dll doit se trouver a cote de
    // l'executable C# (voir GeoTechRadier.csproj, copie post-build).
    private const string DllName = "GeoInfraRadierCalcul.dll";

    [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
    public static extern int GIR_ContourRectangle(double lx, double ly,
        [Out] double[] xyOut, out int nbSommetsOut);

    [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
    public static extern int GIR_ContourL(double lx, double ly, double bx, double by,
        [Out] double[] xyOut, out int nbSommetsOut);

    [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
    public static extern int GIR_ContourU(double lx, double ly, double x1, double x2, double by,
        [Out] double[] xyOut, out int nbSommetsOut);

    [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
    public static extern int GIR_Dimensionner(ref GIR_Entrees entrees, out GIR_Dimensions dims);

    [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
    public static extern int GIR_ObtenirMaillage(ref GIR_Entrees entrees, ref GIR_Maillage maillage);

    [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
    public static extern int GIR_Calculer(ref GIR_Entrees entrees, ref GIR_Resultats resultats);

    [DllImport(DllName, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern int GIR_Exporter(ref GIR_Entrees entrees, string? cheminXlsx,
        string? cheminSvg, string? cheminHtml);

    [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
    public static extern int GIR_DernierMessageErreur(
        [MarshalAs(UnmanagedType.LPArray)] byte[] tampon, int tailleTampon);
}
