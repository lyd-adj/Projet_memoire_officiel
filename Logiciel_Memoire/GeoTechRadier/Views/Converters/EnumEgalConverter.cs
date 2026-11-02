using System.Globalization;
using System.Windows.Data;

namespace GeoTechRadier.Views.Converters;

/// <summary>Lie un RadioButton.IsChecked a une valeur d'enum precise (ex. Forme == FormeRadier.L),
/// pour piloter la selection de la geometrie du radier (section 6 du cahier des charges).</summary>
public class EnumEgalConverter : IValueConverter
{
    public object Convert(object? value, Type targetType, object? parameter, CultureInfo culture)
        => value?.ToString() == parameter?.ToString();

    public object? ConvertBack(object value, Type targetType, object? parameter, CultureInfo culture)
        => (value is bool b && b) ? Enum.Parse(targetType, parameter!.ToString()!) : Binding.DoNothing;
}

/// <summary>Affiche un element uniquement si la forme selectionnee correspond au parametre.</summary>
public class EnumVersVisibiliteConverter : IValueConverter
{
    public object Convert(object? value, Type targetType, object? parameter, CultureInfo culture)
        => value?.ToString() == parameter?.ToString()
            ? System.Windows.Visibility.Visible
            : System.Windows.Visibility.Collapsed;

    public object ConvertBack(object value, Type targetType, object? parameter, CultureInfo culture)
        => throw new NotSupportedException();
}
