#ifndef MY_IMAGE_H
#define MY_IMAGE_H



class my_image {
public:
	const static int bytes_per_pixel = 3;
	my_image()
		: data(nullptr), width(0), height(0), bytes_per_scanline(0) {}
	my_image(const char* filename);

	int get_width() {
		return this->width;
	}
	int get_height() {
		return this->height;
	}
	int get_bytes_per_scanline() {
		return this->bytes_per_scanline;
	}

private:
	unsigned char* data;
	int width, height;
	int bytes_per_scanline;

};

#endif // ! MY_IMAGE_H
