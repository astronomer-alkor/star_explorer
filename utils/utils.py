from math import log10, pi, sqrt


def delta_m(e1, e2):
    return log10(e1 / e2) / 0.4


def get_luminosity(absolute_magnitude):
    return 2.512 ** (4.83 - absolute_magnitude)


def get_mass(luminosity):
    return luminosity ** (1 / 3.9)


def get_temperature(color_index):
    return 7920 / (color_index + 0.72)


def get_age_by_mass(mass):
    return 1E10 / mass ** 3


def get_density(mass, radius):
    return mass / (4 / 3 * pi * radius ** 3)


def get_radius(luminosity, temperature):
    stefan_boltzmann_constant = 5.67036713E-8
    luminosity_of_sun = 3.827E26
    sun_radius = 6.9551E8
    return sqrt(luminosity * luminosity_of_sun / (4 * pi * stefan_boltzmann_constant *
                                                  pow(temperature, 4))) / sun_radius


def get_absolute_magnitude_from_hertzsprung_russell_diagram(color_index):
    color_index *= 10
    if color_index > 19:
        return ((0.1581608 * pow(color_index, 2) + 4.3310835 * color_index + 6.2940268) / 10 +
                ((1.3948753E-6) * pow(color_index, 8.0) - (8.7648505E-5) * pow(color_index, 7.0) +
                 0.001944 * pow(color_index, 6) - 0.0153264 * pow(color_index, 5) - 0.0389616 *
                 pow(color_index, 4) + 1.1514965 * pow(color_index, 3) - 4.7452887 * pow(color_index, 2)
                 + 9.6724702 * color_index + 7.8898473) / 10) / 2

    return ((1.3948753E-6) * pow(color_index, 8.0) - (8.7648505E-5) * pow(color_index, 7.0)
            + 0.001944 * pow(color_index, 6) - 0.0153264 * pow(color_index, 5) - 0.0389616
            * pow(color_index, 4) + 1.1514965 * pow(color_index, 3) - 4.7452887 * pow(color_index, 2)
            + 9.6724702 * color_index + 7.8898473) / 10
