from collections import namedtuple
from typing import List

import astroalign as aa
import matplotlib.pyplot as plt
from astropy.io import fits

from utils.utils import (
    delta_m,
    get_absolute_magnitude_from_hertzsprung_russell_diagram,
    get_luminosity,
    get_mass,
    get_temperature,
    get_age_by_mass,
    get_radius,
    get_density,
)

MIN_BRIGHTNESS_STAR = 2200
ADDITIONAL_BRIGHTNESS = 20000

MIN_PIXELS_ON_STAR = 2
MAX_PIXELS_ON_STAR = 4


def load_image(path):
    with fits.open(path) as hdul:
        image = hdul[0].data
    return image


def adjust_brightness(image, coefficient=1):
    for i, row in enumerate(image):
        for j, pixel in enumerate(row):
            if pixel > MIN_BRIGHTNESS_STAR:
                image[i][j] += ADDITIONAL_BRIGHTNESS * coefficient


Point = namedtuple('Point', ['x', 'y', 'signal'])
Coordinate = namedtuple('Coordinate', ['x', 'y'])


class Star:
    def __init__(self):
        self.brightness = 0
        self.pixel_positions = set()


def get_stars_brightness(image):
    stars = []
    visited = set()
    queue = []
    for i, row in enumerate(image):
        for j, pixel in enumerate(row):
            if (j, i) not in visited:
                visited.add((j, i))
                if pixel >= MIN_BRIGHTNESS_STAR:
                    star = Star()
                    queue.append(Point(j, i, pixel))
                    while queue:
                        tmp = queue.pop(0)
                        star.brightness += tmp.signal
                        star.pixel_positions.add(Coordinate(tmp.x, tmp.y))
                        visited.add((tmp.x, tmp.y))

                        coordinates = ((tmp.y - 1, tmp.x), (tmp.y, tmp.x - 1), (tmp.y, tmp.x + 1), (tmp.y + 1, tmp.x))
                        for y, x in coordinates:
                            if (x, y) not in visited and all((y >= 0, x >= 0, y < len(image), x < len(image[0]))):
                                if image[y][x] >= MIN_BRIGHTNESS_STAR:
                                    visited.add((x, y))
                                    queue.append(Point(x, y, pixel))
                    if len(star.pixel_positions) in range(MIN_PIXELS_ON_STAR, MAX_PIXELS_ON_STAR + 1):
                        stars.append(star)
    return stars


def normalize_arrays(blue_stars: List[Star], visual_stars: List[Star], luminosity_stars: List[Star]):
    normalized_blue_stars = []
    normalized_visual_stars = []
    normalized_luminosity_stars = []
    for blue_star in blue_stars:
        aligned = False
        for visual_star in visual_stars:
            for luminosity_star in luminosity_stars:
                if blue_star.pixel_positions.intersection(visual_star.pixel_positions, luminosity_star.pixel_positions):
                    normalized_blue_stars.append(blue_star)
                    normalized_visual_stars.append(visual_star)
                    normalized_luminosity_stars.append(luminosity_star)
                    aligned = True
                    break
            if aligned:
                break
    return normalized_blue_stars, normalized_visual_stars, normalized_luminosity_stars


class AdjustBrightnessToDetectAlignment:
    def __init__(self, *images):
        self.images = images

    def __enter__(self):
        for image in self.images:
            adjust_brightness(image)

    def __exit__(self, exc_type, exc_val, exc_tb):
        for image in self.images:
            adjust_brightness(image, coefficient=-1)


if __name__ == '__main__':
    visual_image = load_image('images/ngc2808/raw-T30-alkor-ngc2808-20210304-002127-V-BIN1-W-015-001.fit')
    blue_image = load_image('images/ngc2808/raw-T30-alkor-ngc2808-20210304-002043-B-BIN1-W-015-001.fit')
    luminosity_image = load_image('images/ngc2808/raw-T30-alkor-ngc2808-20210304-002208-Luminance-BIN1-W-015-001.fit')
    # known_star_m = 2.5

    with AdjustBrightnessToDetectAlignment(visual_image, blue_image, luminosity_image):
        for image in (blue_image, luminosity_image):
            transformation, _ = aa.find_transform(image, visual_image)
            x, y = map(int, transformation.translation)
            for i, row in enumerate(image):
                for j, item in enumerate(row):
                    try:
                        image[i + y][j + x] = image[i][j]
                    except IndexError:
                        pass

    blue_stars, visual_starts, luminosity_stars = (
        get_stars_brightness(blue_image),
        get_stars_brightness(visual_image),
        get_stars_brightness(luminosity_image)
    )
    print(len(blue_stars), len(visual_starts), len(luminosity_stars))

    blue_stars, visual_stars, luminosity_stars = normalize_arrays(blue_stars, visual_starts, luminosity_stars)
    print(len(blue_stars), len(visual_starts), len(luminosity_stars))

    blue_stars_delta_m = [delta_m(luminosity_star.brightness, blue_star.brightness)
                          for luminosity_star, blue_star in zip(luminosity_stars, blue_stars)]

    visual_stars_delta_m = [delta_m(luminosity_star.brightness, visual_star.brightness)
                            for luminosity_star, visual_star in zip(luminosity_stars, visual_stars)]

    color_indexes = [b - v for b, v in zip(blue_stars_delta_m, visual_stars_delta_m)]
    color_temperatures = [get_temperature(color_index) for color_index in color_indexes]
    absolute_magnitudes = list(map(get_absolute_magnitude_from_hertzsprung_russell_diagram, color_indexes))

    print(color_indexes)
    print(absolute_magnitudes)

    index = int(input('enter index of pivot star: '))

    luminosity = get_luminosity(absolute_magnitudes[index])
    mass = get_mass(luminosity)
    age = get_age_by_mass(mass)
    radius = get_radius(luminosity, color_temperatures[index])
    density = get_density(mass, radius)

    print(f"""
        Mass: {mass}
        Radius: {radius}
        Density: {density}
        Luminosity: {luminosity}
        Age: {age}
    """)
    # adjust_brightness(visual)
    # adjust_brightness(blue)
    # plt.imshow(visual, cmap='gray')
    # plt.show()
    # plt.imshow(blue, cmap='gray')
    # plt.show()
