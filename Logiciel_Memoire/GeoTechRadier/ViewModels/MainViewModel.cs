using System.Collections.ObjectModel;
using System.Linq;
using System.Windows;
using System.Windows.Media;
using CommunityToolkit.Mvvm.ComponentModel;
using CommunityToolkit.Mvvm.Input;
using GeoTechRadier.Models;
using GeoTechRadier.Services;

namespace GeoTechRadier.ViewModels;

/// <summary>Ecrans de l'application, dans l'ordre du workflow decrit au cahier des charges
/// (section 4/35) : Projet -> Radier &amp; Maillage -> Sol -> Charges -> Calcul -> Resultats ->
/// Visualisation. Accueil est l'ecran affiche tant qu'aucun projet n'est ouvert.</summary>
public enum EcranApplication
{
    Accueil,
    Projet,
    RadierMaillage,
    Sol,
    Charges,
    Calcul,
    Resultats,
    Visualisation
}

/// <summary>Orchestre la navigation entre modules et le lancement du calcul. N'effectue
/// elle-meme aucun calcul mecanique : delegue integralement a RadierCalculService, qui appelle
/// le moteur natif GeoInfraRadierCalcul.dll (lui-meme une compilation de GeoInfraRadier.cpp).</summary>
public partial class MainViewModel : ObservableObject
{
    public const string TitreLogiciel = "GEOINFRARADIER";
    public const string SousTitreLogiciel = "Dimensionnement geotechnique des radiers";
    public const string Version = "1.0";

    // Unites reellement utilisees par le moteur (cf. entete de GeoInfraRadier.cpp : "Unites :
    // kN, m => pressions et modules en kPa"). Affichees telles quelles, sans conversion.
    public const string UnitesAffichees = "m - kN - kPa";

    public ProjetViewModel Projet { get; } = new();
    public DonneesViewModel Donnees { get; } = new();
    public ResultatsViewModel Resultats { get; } = new();
    public CalculViewModel Calcul { get; } = new();

    [ObservableProperty] private string messageStatut = "Pret.";
    [ObservableProperty] private bool calculEnCours;
    [ObservableProperty] private bool projetEstOuvert;
    [ObservableProperty] private EcranApplication ecranActif = EcranApplication.Accueil;

    // ---- Ecran d'accueil (section 5) ---------------------------------------
    // Liste honnete : aucun projet recent n'est fabrique. Tant que la persistance de projet
    // (F02/F03) n'est pas implementee, cette liste reste vide plutot que de simuler un
    // historique -- conformement a la consigne de ne jamais creer de fausse fonctionnalite.
    public ObservableCollection<string> ProjetsRecents { get; } = new();

    [RelayCommand]
    private void NouveauProjet()
    {
        ProjetEstOuvert = true;
        EcranActif = EcranApplication.Projet;
        MessageStatut = "Nouveau projet cree.";
    }

    [RelayCommand]
    private void ValiderProjet()
    {
        // Les informations du projet sont deja liees en direct (TwoWay binding) au
        // ProjetViewModel : il n'y a donc rien de plus a "enregistrer" en memoire. Ce bouton
        // confirme explicitement la saisie a l'utilisateur et l'emmene vers l'etape suivante du
        // workflow (section 35 du cahier des charges), plutot que de le laisser deviner qu'il
        // peut simplement changer d'onglet.
        MessageStatut = $"Projet '{Projet.NomProjet}' valide.";
        EcranActif = EcranApplication.RadierMaillage;
    }

    [RelayCommand]
    private void OuvrirProjet()
    {
        // F03 (ouverture d'un projet enregistre) n'est pas encore implementee : le bouton
        // existe (section 7) mais informe clairement l'utilisateur plutot que de simuler une
        // ouverture. Voir section 58 du cahier des charges : ne jamais fabriquer une fausse
        // fonctionnalite.
        MessageBox.Show(
            "L'enregistrement et l'ouverture de projet ne sont pas encore implementes dans cette version.\n" +
            "Utilisez 'Nouveau projet' pour commencer une etude.",
            "Fonctionnalite a venir", MessageBoxButton.OK, MessageBoxImage.Information);
    }

    [RelayCommand]
    private void Naviguer(EcranApplication ecran) => EcranActif = ecran;

    private bool PeutExporter() => Resultats.ADesResultats;

    // Les exports reutilisent directement les fonctions existantes du moteur
    // (exporterIsovaleursTassementSVG / exporterGraphique3DTassementHTML / exporterResultatsExcel,
    // via GIR_Exporter) : la seule grandeur exportable graphiquement aujourd'hui est le
    // tassement w -- aucune autre grandeur (ks, pression, Rz...) n'a de fonction d'export dans
    // GeoInfraRadier.cpp, donc l'IHM ne propose pas de selecteur qui ferait semblant d'en changer.
    [RelayCommand(CanExecute = nameof(PeutExporter))]
    private void GenererCarte2D() => GenererEtOuvrir(svg: true);

    [RelayCommand(CanExecute = nameof(PeutExporter))]
    private void GenererVue3D() => GenererEtOuvrir(svg: false);

    [RelayCommand(CanExecute = nameof(PeutExporter))]
    private void ExporterExcel()
    {
        try
        {
            var dossier = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "GeoInfraRadier");
            System.IO.Directory.CreateDirectory(dossier);
            var chemin = System.IO.Path.Combine(dossier, "resultats_radier.xlsx");
            RadierCalculService.Exporter(Donnees.VersDonneesCalcul(), chemin, null, null);
            System.Diagnostics.Process.Start(new System.Diagnostics.ProcessStartInfo(chemin) { UseShellExecute = true });
            MessageStatut = $"Classeur Excel exporte : {chemin}";
        }
        catch (MoteurCalculException ex)
        {
            MessageBox.Show(ex.Message, "Erreur d'export", MessageBoxButton.OK, MessageBoxImage.Error);
        }
    }

    private void GenererEtOuvrir(bool svg)
    {
        try
        {
            var dossier = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "GeoInfraRadier");
            System.IO.Directory.CreateDirectory(dossier);
            var cheminSvg = System.IO.Path.Combine(dossier, "isovaleurs_tassement.svg");
            var cheminHtml = System.IO.Path.Combine(dossier, "tassement_3D.html");
            RadierCalculService.Exporter(Donnees.VersDonneesCalcul(), null,
                svg ? cheminSvg : null, svg ? null : cheminHtml);
            var chemin = svg ? cheminSvg : cheminHtml;
            System.Diagnostics.Process.Start(new System.Diagnostics.ProcessStartInfo(chemin) { UseShellExecute = true });
            MessageStatut = $"Vue generee et ouverte : {chemin}";
        }
        catch (MoteurCalculException ex)
        {
            MessageBox.Show(ex.Message, "Erreur d'export", MessageBoxButton.OK, MessageBoxImage.Error);
        }
    }

    // Taille du canevas d'apercu (px) dans les vues -- doit rester coherente avec le Canvas
    // Width/Height defini dans le XAML.
    private const double TailleApercuPx = 280.0;

    [RelayCommand]
    private void ActualiserApercu()
    {
        try
        {
            var sommets = RadierCalculService.ConstruireContour(Donnees.VersDonneesCalcul());
            Donnees.ContourApercu.Clear();
            foreach (var s in sommets) Donnees.ContourApercu.Add(s);

            // Mise a l'echelle simple dans le canevas, avec inversion de l'axe Y (repere ecran
            // vers le bas, repere du radier vers le haut) : pure presentation, la geometrie
            // reelle (sommets) reste celle renvoyee par le moteur.
            double maxDim = Math.Max(Donnees.Lx, Donnees.Ly);
            if (maxDim <= 0) maxDim = 1.0;
            Donnees.ContourPoints = new PointCollection(sommets.Select(s =>
                new Point(s.X / maxDim * TailleApercuPx, TailleApercuPx - s.Y / maxDim * TailleApercuPx)));

            Donnees.Maillage = RadierCalculService.Dimensionner(Donnees.VersDonneesCalcul());

            // Grille reelle du maillage (section 15/16) : un quadrilatere par element actif,
            // en coordonnees deja mises a l'echelle (memes conventions que le contour ci-dessus).
            var geometrie = RadierCalculService.ObtenirMaillageGeometrie(Donnees.VersDonneesCalcul());
            Donnees.ElementsAffichage.Clear();
            foreach (var quad in geometrie.Elements)
            {
                Donnees.ElementsAffichage.Add(new PointCollection(quad.Select(s =>
                    new Point(s.X / maxDim * TailleApercuPx, TailleApercuPx - s.Y / maxDim * TailleApercuPx))));
            }

            const double rayon = 5.0;
            Donnees.PointsAffichage.Clear();
            foreach (var p in Donnees.Ponctuelles)
            {
                Donnees.PointsAffichage.Add(new Point(
                    p.X / maxDim * TailleApercuPx - rayon,
                    TailleApercuPx - p.Y / maxDim * TailleApercuPx - rayon));
            }
            MessageStatut = $"Contour a jour : {sommets.Count} sommets, {Donnees.Maillage.NbElements} elements, " +
                $"{Donnees.Maillage.Nn} noeuds.";
        }
        catch (MoteurCalculException ex)
        {
            MessageStatut = $"Geometrie invalide : {ex.Message}";
        }
    }

    private bool PeutCalculer() => !CalculEnCours;

    [RelayCommand(CanExecute = nameof(PeutCalculer))]
    private async Task LancerCalculAsync()
    {
        var erreurs = Donnees.Valider();
        if (erreurs.Count > 0)
        {
            MessageBox.Show(
                "Impossible de lancer le calcul :\n\n- " + string.Join("\n- ", erreurs),
                "Donnees incompletes ou incoherentes", MessageBoxButton.OK, MessageBoxImage.Warning);
            return;
        }

        CalculEnCours = true;
        EcranActif = EcranApplication.Calcul;
        Calcul.Reinitialiser();
        Calcul.EnCours = true;
        MessageStatut = "Calcul en cours...";
        Resultats.Vider();
        try
        {
            var donneesCalcul = Donnees.VersDonneesCalcul();
            // Le calcul (assemblage + resolution du systeme) peut etre couteux sur un maillage
            // fin : execute hors du thread UI (section 37.3) pour ne jamais geler l'interface.
            var resultats = await Task.Run(() => RadierCalculService.Calculer(donneesCalcul, etape =>
            {
                // Callback invoque depuis le thread de calcul : re-synchronise sur le thread UI.
                Application.Current.Dispatcher.Invoke(() =>
                {
                    int index = etape == RadierCalculService.EtapeNative.PreparationMaillage ? 0 : 1;
                    if (index == 1) Calcul.MarquerEtape(0, EtatEtape.Termine);
                    Calcul.MarquerEtape(index, EtatEtape.EnCours);
                });
            }));

            Calcul.MarquerEtape(1, EtatEtape.Termine);
            Calcul.EnCours = false;
            Calcul.TermineAvecSucces = true;

            Resultats.Charger(resultats);
            MessageStatut = $"Calcul termine : {resultats.Nn} noeuds, {resultats.NbElements} elements actifs. " +
                $"Equilibre : somme reactions = {resultats.SommeReactionP:0.###} kN, somme Rz = {resultats.SommeRz:0.###} kN.";
            EcranActif = EcranApplication.Resultats;
        }
        catch (MoteurCalculException ex)
        {
            Calcul.EnCours = false;
            Calcul.MessageErreur = ex.Message;
            foreach (var e in Calcul.Etapes) if (e.Etat == EtatEtape.EnCours) e.Etat = EtatEtape.Erreur;
            MessageStatut = "Erreur du moteur de calcul.";
            MessageBox.Show(ex.Message, "Erreur du moteur de calcul", MessageBoxButton.OK, MessageBoxImage.Error);
        }
        catch (Exception ex)
        {
            Calcul.EnCours = false;
            Calcul.MessageErreur = ex.Message;
            foreach (var e in Calcul.Etapes) if (e.Etat == EtatEtape.EnCours) e.Etat = EtatEtape.Erreur;
            MessageStatut = "Erreur inattendue.";
            MessageBox.Show(ex.Message, "Erreur inattendue", MessageBoxButton.OK, MessageBoxImage.Error);
        }
        finally
        {
            CalculEnCours = false;
        }
    }
}
