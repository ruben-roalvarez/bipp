#ifndef PROFILE_WORLOADS_HPP
#define PROFILE_WORLOADS_HPP

#include <fstream>
#include <iostream>
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>

static int workload_number = 0;

/**
 * Finds even integer not less than n, with prime factors no larger than 5 (ie, "smooth").
 * Doesn't need to be done in HW because the FFT has to be a power of 2 (by definition it only has mutiples of 2).
 * If used, only apply it to numbers less than 1e12.
 */
int next235(int n)
{
  if (n <= 2) return 2;
  if (n % 2 == 1) n += 1;                   // even
  int nplus  = n - 2;                       // to cancel out the +=2 at start of loop
  int numdiv = 2;                           // a dummy that is >1
  while (numdiv > 1) {
    nplus += 2;                             // stays even
    numdiv = nplus;
    while (numdiv % 2 == 0) numdiv /= 2;    // remove all factors of 2,3,5...
    while (numdiv % 3 == 0) numdiv /= 3;
    while (numdiv % 5 == 0) numdiv /= 5;
  }
  return nplus;
}

template <typename T>
double calculate_width_out(T *points, int size) {
    double high = (double)*std::max_element(points, points + size);
    double low = (double)*std::min_element(points, points + size);
    double width_out = (high - low) / 2.0;
    double center_out = (high + low) / 2.0;
    if (std::abs(center_out) < 0.1 * (width_out)) {
        width_out += std::abs(center_out);
        center_out = 0.0;
    }
    return width_out;
}

template <typename T>
int get_fft_size(T *pixel_points, int numPixels, T *uvw_points, int numFreqs, int nspread){
    
    constexpr double PI = M_PI;

    double width_out = calculate_width_out<T>(pixel_points, numPixels);
    double width_in = calculate_width_out<T>(uvw_points, numFreqs);

    // ensures width_out * width_in >= 1
    double XS = width_out * width_in;
    if (XS < 1.0) XS = 1.0;

    // Initial value
    int grid_size = 4.0 * XS / PI + nspread + 1;
    if (grid_size < 2 * nspread) grid_size = 2 * nspread; // Minimum
    grid_size = next235(grid_size)*2;

    return grid_size;
}

template <typename T>
void read_from_file(const std::string& filename, std::vector<std::complex<T>> &complex_vector) {
  std::ifstream file(filename, std::ios::binary);
  
  if (file.is_open()) {
    T real, imag;
    while (file.read(reinterpret_cast<char*>(&real), sizeof(T)) &&
           file.read(reinterpret_cast<char*>(&imag), sizeof(T))) {
      complex_vector.emplace_back(real, imag);
    }
    file.close();
  }else {
    std::cerr << "Unable to open file" << std::endl;
  }
  
  return;
}

template <typename T>
void read_from_file(const std::string& filename, std::vector<T> &vector) {
  std::ifstream file(filename, std::ios::binary);
  
  if (file.is_open()) {
    T value;
    while (file.read(reinterpret_cast<char*>(&value), sizeof(T))) {
      vector.emplace_back(value);
    }
    file.close();
  }else {
    std::cerr << "Unable to open file" << std::endl;
  }
  
  return;
}

template <typename T>
void write_to_file(const std::string& filename, std::complex<T> * complex_vector, size_t size) {
  std::ofstream file(filename, std::ios::binary);

  if (file.is_open()) {
    for (size_t i = 0; i < size; i++) {
      T real = std::real(complex_vector[i]);
      T imag = std::imag(complex_vector[i]);
      file.write(reinterpret_cast<const char*>(&real), sizeof(T));
      file.write(reinterpret_cast<const char*>(&imag), sizeof(T));
    }
    file.close();
  }
}

template <typename T>
void write_to_file(const std::string& filename, T * vector, size_t size) {
  std::ofstream file(filename, std::ios::binary);

  if (file.is_open()) {
    for (size_t i = 0; i < size; i++) {
      T real = vector[i];
      file.write(reinterpret_cast<const char*>(&real), sizeof(T));
    }
    file.close();
  }
}

#endif // PROFILE_WORLOADS_HPP
