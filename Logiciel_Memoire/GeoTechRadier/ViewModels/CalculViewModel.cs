using System.Collections.ObjectModel;
using System.Linq;
using CommunityToolkit.Mvvm.ComponentModel;

namespace GeoTechRadier.ViewModels;

public enum EtatEtape { EnAttente, EnCours, Termine, Erreur }

public partial class EtapeCalcul : ObservableObject
{
    [ObservableProperty] private string libelle = string.Empty;
    [ObservableProperty] private EtatEtape etat = EtatEtape.EnAttente;
}

/// <summary>Suivi du module "Calcul" (section 14/16 du cahier des charges). Le contrat natif
/// actuel (GIR_Dimensionner puis GIR_Calculer) n'expose que DEUX etapes reellement observables
/// depuis l'exterieur du moteur -- la finesse en 8 etapes decrite dans le cahier des charges
/// (Ke, Kplaque, Ts, Mp/B, KT, FP, resolution, post-traitement) correspond a des operations
/// internes a solveSsi() qui ne remontent aujourd'hui aucun evenement de progression.
/// Afficher ces 8 etapes ferait croire a une progression que l'IHM n'observe pas reellement :
/// on affiche donc honnetement les deux etapes reelles, en gardant la porte ouverte pour un
/// futur enrichissement du contrat GIR_* si une instrumentation plus fine est ajoutee au moteur.</summary>
public partial class CalculViewModel : ObservableObject
{
    public ObservableCollection<EtapeCalcul> Etapes { get; } = new()
    {
        new EtapeCalcul { Libelle = "Preparation du maillage (dimensionnement du systeme)" },
        new EtapeCalcul { Libelle = "Assemblage et resolution du systeme couple radier-sol" },
    };

    [ObservableProperty] private bool enCours;
    [ObservableProperty] private bool termineAvecSucces;
    [ObservableProperty] private string? messageErreur;
    [ObservableProperty] private double progressionPourcent;

    public void Reinitialiser()
    {
        foreach (var e in Etapes) e.Etat = EtatEtape.EnAttente;
        ProgressionPourcent = 0;
        MessageErreur = null;
        TermineAvecSucces = false;
    }

    public void MarquerEtape(int index, EtatEtape etat)
    {
        if (index < 0 || index >= Etapes.Count) return;
        Etapes[index].Etat = etat;
        var terminees = Etapes.Count(e => e.Etat == EtatEtape.Termine);
        ProgressionPourcent = 100.0 * terminees / Etapes.Count;
    }
}
