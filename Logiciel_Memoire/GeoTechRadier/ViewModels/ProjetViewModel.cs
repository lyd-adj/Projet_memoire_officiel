using CommunityToolkit.Mvvm.ComponentModel;

namespace GeoTechRadier.ViewModels;

/// <summary>Informations generales du projet (section 16 du cahier des charges). Volontairement
/// minimal tant que le calcul n'est pas integre et valide -- pas de sauvegarde/chargement de
/// projet a ce stade (Lot 3).</summary>
public partial class ProjetViewModel : ObservableObject
{
    [ObservableProperty] private string nomProjet = "Nouveau projet";
    [ObservableProperty] private string reference = string.Empty;
    [ObservableProperty] private string description = string.Empty;
    [ObservableProperty] private string ingenieur = string.Empty;
    [ObservableProperty] private DateTime date = DateTime.Now;
}
