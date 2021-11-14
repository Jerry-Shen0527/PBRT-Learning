#include "texture.h"
#include "rtw_stb_image.h"


image_texture::image_texture(const char* filename) {
    auto components_per_pixel = bytes_per_pixel;

    data = stbi_load(
        filename, &width, &height, &components_per_pixel, components_per_pixel);

    if (!data) {
        std::cerr << "ERROR: Could not load texture image file '" << filename << "'.\n";
        width = height = 0;
    }
    else
    {
        //std::cout << "load image file successfully!\n" << "width: " << width << "    height:" << height << std::endl;
    }
    bytes_per_scanline = bytes_per_pixel * width;
}