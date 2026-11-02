using System.Collections.ObjectModel;
using System.Windows.Media;
using CommunityToolkit.Mvvm.ComponentModel;
using CommunityToolkit.Mvvm.Input;
using GeoTechRadier.Models;
using GeoTechRadier.Services;

namespace GeoTechRadier.ViewModels;

/// <summary>Saisie complete d'un cas de calcul : geometrie, maillage, materiau, sol, chargements.
/// Miroir exact, ecran par ecran, de la sequence de saisie de main() dans GeoInfraRadier.cpp
/// (section 6 a 15 du cahier des charges) -- aucun parametre n'est invente ici.</summary>
public partial class DonneesViewModel : ObservableObject
{
    // ---- 5.1 / 6-9 : Geometrie -------------------------------------------------
    [ObservableProperty] private FormeRadier forme = FormeRadier.Rectangulaire;

    [ObservableProperty] private double lx = 5.0;
    [ObservableProperty] private double ly = 5.0;

    [ObservableProperty] private double bx = 2.5;   // echancrure L
    [ObservableProperty] private double by = 2.5;   // echancrure L

    [ObservableProperty] private double x1 = 1.25;  // echancrure U
    [ObservableProperty] private double x2 = 3.75;  // echancrure U
    [ObservableProperty] private double byU = 2.5;  // echancrure U

    [ObservableProperty] private double tailleX = 0.5;
    [ObservableProperty] private double tailleY = 0.5;

    // Apercu du contour (rafraichi a la demande, cf. MainViewModel.ActualiserApercu) : utilise
    // pour la visualisation avant calcul (section 10) sans dupliquer la geometrie native --
    // les sommets viennent de RadierCalculService.ConstruireContour (appel direct au moteur).
    public ObservableCollection<Sommet> ContourApercu { get; } = new();

    // Points en pixels (deja mis a l'echelle) pour le dessin du Polygon dans DonneesView --
    // recalcules par MainViewModel.ActualiserApercu a partir de ContourApercu.
    [ObservableProperty] private PointCollection contourPoints = new();

    // Indicateurs de maillage (section 8/16), rafraichis par MainViewModel.ActualiserApercu via
    // RadierCalculService.Dimensionner -- disponibles avant meme de lancer un calcul complet.
    [ObservableProperty] private MaillageInfo maillage = new();

    // Contour de chaque element du maillage (deja mis a l'echelle, en pixels), pour dessiner la
    // grille reelle (section 15/16) -- rempli par MainViewModel.ActualiserApercu a partir de
    // RadierCalculService.ObtenirMaillageGeometrie (donnees reelles du moteur, pas une grille
    // approximative recalculee cote C#).
    public ObservableCollection<PointCollection> ElementsAffichage { get; } = new();

    // Positions (deja mises a l'echelle, en pixels) des charges ponctuelles dans l'apercu
    // graphique du module Charges -- recalculees par MainViewModel.ActualiserApercu, a partir
    // des memes coordonnees reelles (Ponctuelles) que celles transmises au moteur.
    public ObservableCollection<System.Windows.Point> PointsAffichage { get; } = new();

    // ---- 12 : Materiau du radier ------------------------------------------------
    [ObservableProperty] private double eBeton = 3.1e7;   // kPa
    [ObservableProperty] private double nuBeton = 0.25;
    [ObservableProperty] private double hRadier = 0.4;    // m

    // ---- 13 : Sol multicouche ----------------------------------------------------
    public ObservableCollection<CoucheSol> Couches { get; } = new() { new CoucheSol() };
    [ObservableProperty] private CoucheSol? coucheSelectionnee;

    [RelayCommand]
    private void AjouterCouche() => Couches.Add(new CoucheSol());

    [RelayCommand]
    private void SupprimerCouche()
    {
        if (CoucheSelectionnee != null) Couches.Remove(CoucheSelectionnee);
    }

    // ---- 14 : Chargements ---------------------------------------------------------
    public ObservableCollection<ChargeSurfacique> Surfaces { get; } = new();
    [ObservableProperty] private ChargeSurfacique? surfaceSelectionnee;

    [RelayCommand]
    private void AjouterSurface() => Surfaces.Add(new ChargeSurfacique());

    [RelayCommand]
    private void SupprimerSurface()
    {
        if (SurfaceSelectionnee != null) Surfaces.Remove(SurfaceSelectionnee);
    }

    public ObservableCollection<ChargePonctuelle> Ponctuelles { get; } = new()
    {
        new ChargePonctuelle { X = 2.5, Y = 2.5 }
    };
    [ObservableProperty] private ChargePonctuelle? ponctuelleSelectionnee;

    [RelayCommand]
    private void AjouterPonctuelle() => Ponctuelles.Add(new ChargePonctuelle());

    [RelayCommand]
    private void SupprimerPonctuelle()
    {
        if (PonctuelleSelectionnee != null) Ponctuelles.Remove(PonctuelleSelectionnee);
    }

    public ObservableCollection<ChargeLineique> Lineiques { get; } = new();
    [ObservableProperty] private ChargeLineique? lineiqueSelectionnee;

    [RelayCommand]
    private void AjouterLineique() => Lineiques.Add(new ChargeLineique());

    [RelayCommand]
    private void SupprimerLineique()
    {
        if (LineiqueSelectionnee != null) Lineiques.Remove(LineiqueSelectionnee);
    }

    /// <summary>Convertit la saisie en DTO consomme par RadierCalculService.</summary>
    public DonneesCalcul VersDonneesCalcul()
    {
        var d = new DonneesCalcul
        {
            Forme = Forme,
            Lx = Lx,
            Ly = Ly,
            Bx = Bx,
            By = By,
            X1 = X1,
            X2 = X2,
            ByU = ByU,
            TailleX = TailleX,
            TailleY = TailleY,
            EBeton = EBeton,
            NuBeton = NuBeton,
            HRadier = HRadier,
        };
        d.Couches.Clear(); d.Couches.AddRange(Couches);
        d.Surfaces.Clear(); d.Surfaces.AddRange(Surfaces);
        d.Ponctuelles.Clear(); d.Ponctuelles.AddRange(Ponctuelles);
        d.Lineiques.Clear(); d.Lineiques.AddRange(Lineiques);
        return d;
    }

    /// <summary>Validation cote IHM, miroir des controles deja presents dans GeoInfraRadier.cpp
    /// (section 7 du cahier des charges) : bloquer AVANT d'appeler le moteur plutot que de
    /// laisser une entree incoherente produire un maillage vide ou un resultat non pertinent.</summary>
    public IReadOnlyList<string> Valider()
    {
        var erreurs = new List<string>();
        if (Lx <= 0 || Ly <= 0) erreurs.Add("Les dimensions Lx et Ly doivent etre strictement positives.");
        if (TailleX <= 0 || TailleY <= 0) erreurs.Add("La taille d'element doit etre strictement positive.");
        if (Forme == FormeRadier.L && !(Bx > 0 && Bx < Lx && By > 0 && By < Ly))
            erreurs.Add("Echancrure en L invalide : il faut 0 < Bx < Lx et 0 < By < Ly.");
        if (Forme == FormeRadier.U && !(X1 >= 0 && X1 < X2 && X2 <= Lx && ByU > 0 && ByU < Ly))
            erreurs.Add("Echancrure en U invalide : il faut 0 <= x1 < x2 <= Lx et 0 < By < Ly.");
        if (HRadier <= 0) erreurs.Add("L'epaisseur du radier doit etre strictement positive.");
        if (EBeton <= 0) erreurs.Add("Le module de Young du beton doit etre strictement positif.");
        if (Couches.Count == 0) erreurs.Add("Au moins une couche de sol est requise.");
        foreach (var c in Couches)
            if (c.E <= 0) erreurs.Add("Le module de Young d'une couche de sol doit etre strictement positif.");
        foreach (var s in Surfaces)
            if (s.Xmax <= s.Xmin || s.Ymax <= s.Ymin)
                erreurs.Add("Une charge surfacique a des bornes invalides (xmax > xmin et ymax > ymin requis).");
        return erreurs;
    }
}
