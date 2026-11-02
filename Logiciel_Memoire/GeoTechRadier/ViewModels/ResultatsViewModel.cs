using System.Collections.ObjectModel;
using CommunityToolkit.Mvvm.ComponentModel;
using GeoTechRadier.Models;

namespace GeoTechRadier.ViewModels;

/// <summary>Presente les resultats renvoyes par le moteur (section 21/22 du cahier des charges).
/// Ne recalcule rien : se contente d'exposer, sous forme de collections bindables aux DataGrid,
/// exactement ce que RadierCalculService.Calculer a rapporte.</summary>
public partial class ResultatsViewModel : ObservableObject
{
    [ObservableProperty] private bool aDesResultats;
    [ObservableProperty] private int nn;
    [ObservableProperty] private int ndof;
    [ObservableProperty] private int nbElements;
    [ObservableProperty] private int nx;
    [ObservableProperty] private int ny;
    [ObservableProperty] private double sommeReactionP;
    [ObservableProperty] private double sommeRz;

    public ObservableCollection<NoeudResultat> Noeuds { get; } = new();
    public ObservableCollection<ElementResultat> Elements { get; } = new();

    [ObservableProperty] private ValeursExtremes extremes = new();

    public void Charger(ResultatsCalcul r)
    {
        Nn = r.Nn; Ndof = r.Ndof; NbElements = r.NbElements; Nx = r.Nx; Ny = r.Ny;
        SommeReactionP = r.SommeReactionP;
        SommeRz = r.SommeRz;

        Noeuds.Clear();
        foreach (var n in r.Noeuds) Noeuds.Add(n);
        Elements.Clear();
        foreach (var e in r.Elements) Elements.Add(e);

        Extremes = ValeursExtremes.Calculer(r);
        ADesResultats = true;
    }

    public void Vider()
    {
        Noeuds.Clear();
        Elements.Clear();
        ADesResultats = false;
    }
}
