using System.Globalization;
using System.Windows.Data;

namespace GeoTechRadier.Views.Converters;

/// <summary>Convertit une epaisseur de couche (m) en hauteur de dessin (px) pour la coupe
/// stratigraphique verticale du module Sol -- pure presentation, la donnee reelle (epaisseur en
/// m) reste celle saisie/transmise au moteur, cette hauteur ne sert qu'au schema.</summary>
public class EpaisseurVersHauteurConverter : IValueConverter
{
    private const double PixelsParMetre = 14.0;
    private const double HauteurMin = 26.0;
    private const double HauteurMax = 140.0;

    public object Convert(object? value, Type targetType, object? parameter, CultureInfo culture)
    {
        double epaisseur = value is double d ? d : 0.0;
        if (epaisseur <= 0) return HauteurMin;
        return Math.Clamp(epaisseur * PixelsParMetre, HauteurMin, HauteurMax);
    }

    public object ConvertBack(object value, Type targetType, object? parameter, CultureInfo culture)
        => throw new NotSupportedException();
}
