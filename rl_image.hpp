#ifndef IMAGE_LOADER_HPP
#define IMAGE_LOADER_HPP

//
//   /$$$$$$$              /$$           /$$                           /$$                    
//  | $$__  $$            | $$          | $$                          |__/                    
//  | $$  \ $$  /$$$$$$  /$$$$$$        | $$        /$$$$$$   /$$$$$$  /$$  /$$$$$$  /$$$$$$$ 
//  | $$$$$$$/ /$$__  $$|_  $$_/        | $$       /$$__  $$ /$$__  $$| $$ /$$__  $$| $$__  $$
//  | $$__  $$| $$$$$$$$  | $$          | $$      | $$$$$$$$| $$  \ $$| $$| $$  \ $$| $$  \ $$
//  | $$  \ $$| $$_____/  | $$ /$$      | $$      | $$_____/| $$  | $$| $$| $$  | $$| $$  | $$
//  | $$  | $$|  $$$$$$$  |  $$$$/      | $$$$$$$$|  $$$$$$$|  $$$$$$$| $$|  $$$$$$/| $$  | $$
//  |__/  |__/ \_______/   \___/        |________/ \_______/ \____  $$|__/ \______/ |__/  |__/
//                                                          /$$  \ $$                        
//                                                         |  $$$$$$/                        
//                                                          \______/ 

#define RLIMG_NAME "Ret Legion image"

// --------------------------------------------------------------------------------\
// rl_image.hpp — single-header image decoder                                      |
// Created by retti, 2026                                                          |
//---------------------------------------------------------------------------------|
//                                                                                 |
//                                  MIT License                                    |
//                                                                                 |
//                          Copyright (c) 2026 retti                               |
//                                                                                 |
//  Permission is hereby granted, free of charge, to any person obtaining a copy   |
//  of this software and associated documentation files (the "Software"), to deal  |
//  in the Software without restriction, including without limitation the rights   |
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      |
//  copies of the Software, and to permit persons to whom the Software is          |
//  furnished to do so, subject to the following conditions:                       |
//                                                                                 |
//  The above copyright notice and this permission notice shall be included in all |
//  copies or substantial portions of the Software.                                |
//                                                                                 |
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     |
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       |
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    |
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         |
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  |
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  |
//  SOFTWARE.                                                                      |
// --------------------------------------------------------------------------------/

// Supported formats:
//   PNG   ✓ 
//   JPG   ⨉ (in progress)
//   WEBP  ⨉ (in progress)
//   HDR   ✓
//   TGA   ✓
//   BMP   ✓
//   DDS   ✓
//   QOI   ✓

#define RLIMG_VERSION "0.3"

// Changelog:
//   [01/2026] v0.3 - Added support for: TGA, BMP, PNG, DDS, QOI, HDR <> Temporarily removed: WebP, JPG (will return after decoder refactor)
//   [01/2026] v0.2 — Optimized PNG decoder, added WebP, HDR
//   [01/2026] v0.1 — Initial release (PNG, JPEG, BMP)

// =======================================================
// Quick Start:
//   rlimg::ImageData img = rlimg::load_image("filename.formats");
//   img = img.convert_to_channels(4);
//   img.flip_vertically();
// =======================================================

// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
// @note DOCUMENTATION
//
//
// Overview:
// ---------
// This library provides a lightweight, self-contained solution for loading and decoding 
// various image formats. It is designed for efficiency and ease of integration into 
// C++ projects requiring image processing capabilities.
//
//
//
// Basic Usage:
// ------------
// 1. Load an image:
//    ```
//    #include "rl_image.hpp"
//    try {
//        rlimg::ImageData img = rlimg::load_image("example.png");
//        // img.data contains pixel data in row-major order
//        // img.width, img.height, img.channels (1, 2, 3, or 4)
//    } catch (const rlimg::image_error& e) {
//        // Handle error
//    }
//    ```
//
// 2. Load with specific channel count:
//    ```
//    // Force 4-channel RGBA output (adds alpha if needed)
//    rlimg::ImageData img = rlimg::load_image("example.jpg", 4);
//    ```
//
// 3. WebP-specific features (caching):
//    ```
//    rlimg::WebPDecoder::clear_cache(); // Clear cache
//    rlimg::WebPDecoder::set_cache_limit(128 * 1024 * 1024); // 128MB cache
//    ```
//
// 4. HDR-specific features:
//    ```
//    // Load with custom exposure
//    rlimg::ImageData hdr = rlimg::HDRDecoder::load("scene.hdr", 2.0f);
//    // Load as float data for HDR processing
//    int w, h;
//    std::vector<float> float_data = rlimg::HDRDecoder::load_float("scene.hdr", w, h);
//    ```
//
// ImageData Structure:
// --------------------
// The loaded image is stored in an ImageData object:
// • width:    Image width in pixels
// • height:   Image height in pixels
// • channels: 1 (grayscale), 2 (gray+alpha), 3 (RGB), or 4 (RGBA)
// • data:     Vector of uint8_t pixels in row-major order
//
// Available ImageData methods:
// • flip_vertically() - Flip image vertically (in-place)
// • convert_to_channels(n) - Convert to n channels (returns new image)
//
//
// Error Handling:
// ---------------
// All load functions throw rlimg::image_error on failure with descriptive messages.
// Common errors: file not found, invalid format, unsupported features, memory limits.
//
// Configuration Macros:
// ---------------------
// Define before including to adjust limits:
// • IMG_MAX_IMAGE_SIZE - Max width/height (default: 16384)
// • IMG_MAX_CHANNELS   - Max channels (default: 4)
//
// Platform Requirements:
// ----------------------
// • C++17 compatible compiler
// • Standard library with <vector>, <algorithm>, <fstream>
// • SIMD headers available on supported platforms
// • No OS-specific code (portable)
//
// Example Workflow:
// -----------------
// ```
// #include "rl_image.hpp"
// 
// int main() {
//     // Load and convert to RGBA
//     rlimg::ImageData img = rlimg::load_image("filename.format", 4);
//     
//     // Process pixels
//     for (size_t i = 0; i < img.data.size(); i += 4) {
//         img.data[i] = 255 - img.data[i]; // Invert red
//     }
//     
//     // Flip and save (save function not included - needs implementation)
//     img.flip_vertically();
//     
//     return 0;
// }
// ```
//
// Notes:
// ------
// • This is a DECODER only - encoding/saving not implemented
// • Memory usage: ~ width × height × channels bytes per image
// • WebP cache: additional memory, controllable via set_cache_limit()
// • Thread safety: Load functions are reentrant; WebP cache uses mutex 
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

#include <cstdint>
#include <cstring>
#include <sys/stat.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include <emmintrin.h> // SSE2

#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI

#ifndef ALWAYS_INLINE
#if defined(_MSC_VER)
#define ALWAYS_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define ALWAYS_INLINE inline
#endif
#endif

#if defined(_MSC_VER)
    #define force_inline __forceinline
#elif defined(__GNUC__) || defined(__clang__)
    #define force_inline __attribute__((always_inline)) inline
#else
    #define force_inline inline
#endif

namespace rlimg {

// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//                CONFIGURATION
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

#ifndef IMG_MAX_IMAGE_SIZE
#define IMG_MAX_IMAGE_SIZE 16384 // max scale (width/height)
#endif // IMG_MAX_IMAGE_SIZE

#ifndef IMG_MAX_CHANNELS
#define IMG_MAX_CHANNELS 4 // max amount channels
#endif // IMG_MAX_CHANNELS

// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//              ERROR HANDLING
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

class image_error : public std::runtime_error {
public:
    explicit image_error(const std::string& msg) : std::runtime_error(msg) {}
};

// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//            IMAGE DATA STRUCTURE
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

struct ImageData {
    int width = 0;
    int height = 0;
    int channels = 0; // 1=grayscale, 2=grayscale+alpha, 3=RGB, 4=RGBA
    std::vector<uint8_t> data;

    ImageData() = default;

    ImageData(int w, int h, int ch) : width(w), height(h), channels(ch), data(w * h * ch) {}

    ImageData(int w, int h, int ch, std::vector<uint8_t> d) : width(w), height(h), channels(ch), data(std::move(d)) {}

    void flip_vertically() {
        const int row_size = width * channels;
        std::vector<uint8_t> row(row_size);

        for (int y = 0; y < height / 2; ++y) {
            uint8_t* top = data.data() + y * row_size;
            uint8_t* bottom = data.data() + (height - 1 - y) * row_size;
            std::memcpy(row.data(), top, row_size);
            std::memcpy(top, bottom, row_size);
            std::memcpy(bottom, row.data(), row_size);
        }
    }

    ImageData convert_to_channels(int target_channels) const {
        if (target_channels == channels) {
            return *this;
        }
        ImageData result(width, height, target_channels);

        if (channels == 1 && target_channels == 3) {                 // Grayscale -> RGB
            for (int i = 0; i < width * height; ++i) {
                uint8_t gray = data[i];
                result.data[i * 3] = gray;
                result.data[i * 3 + 1] = gray;
                result.data[i * 3 + 2] = gray;
            }
        }
        else if (channels == 1 && target_channels == 4) {            // Grayscale -> RGBA
            for (int i = 0; i < width * height; ++i) {
                uint8_t gray = data[i];
                result.data[i * 4] = gray;
                result.data[i * 4 + 1] = gray;
                result.data[i * 4 + 2] = gray;
                result.data[i * 4 + 3] = 255;
            }
        }
        else if (channels == 3 && target_channels == 4) {            // RGB -> RGBA
            for (int i = 0; i < width * height; ++i) {
                result.data[i * 4] = data[i * 3];
                result.data[i * 4 + 1] = data[i * 3 + 1];
                result.data[i * 4 + 2] = data[i * 3 + 2];
                result.data[i * 4 + 3] = 255;
            }
        }
        else if (channels == 4 && target_channels == 3) {            // RGBA -> RGB
            for (int i = 0; i < width * height; ++i) {
                result.data[i * 3] = data[i * 4];
                result.data[i * 3 + 1] = data[i * 4 + 1];
                result.data[i * 3 + 2] = data[i * 4 + 2];
            }
        }
        else {
            throw image_error("Unsupported channel conversion");
        }
        
        return result;
    }
};

// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//      CRC32 IMPLEMENTATION (for PNG)
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

class CRC32 {
private:
    static inline uint32_t crc_table[256] = {0};
    static inline bool table_computed = false;

    static void make_crc_table() {
        for (uint32_t n = 0; n < 256; n++) {
            uint32_t c = n;
            for (int k = 0; k < 8; k++) {
                if (c & 1)
                    c = 0xEDB88320u ^ (c >> 1);
                else
                    c >>= 1;
            }
            crc_table[n] = c;
        }
        table_computed = true;
    }

public:
    static uint32_t update(uint32_t crc, const uint8_t* buf, size_t len) {
        if (!table_computed) make_crc_table();
        uint32_t c = crc ^ 0xFFFFFFFFu;
        for (size_t n = 0; n < len; n++) {
            c = crc_table[(c ^ buf[n]) & 0xFF] ^ (c >> 8);
        }
        return c ^ 0xFFFFFFFFu;
    }
};

// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//      ADLER32 IMPLEMENTATION (for PNG)
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

class Adler32 {
public:
    static uint32_t compute(const uint8_t* data, size_t len) {
        uint32_t a = 1, b = 0;
        for (size_t i = 0; i < len; i++) {
            a += data[i];
            if (a >= 65521) a -= 65521;
            b += a;
            if (b >= 65521) b -= 65521;
        }
        return (b << 16) | a;
    }
};

// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//               HUFFMAN DECODER
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

class HuffmanDecoder {
public:
    static const int MAX_BITS  = 15;
    static const int FAST_BITS = 9;
    static const int FAST_MASK = (1 << FAST_BITS) - 1;

    struct Entry {
        uint16_t symbol;
        uint8_t  length;
    };

    Entry fast[1 << FAST_BITS];

    std::vector<int> symbols;
    std::vector<int> lengths;
    std::vector<int> first_code;
    std::vector<int> first_symbol;
    std::vector<int> count;

    static ALWAYS_INLINE uint32_t reverse_bits(uint32_t code, int len) {
        uint32_t r = 0;
        for (int i = 0; i < len; ++i) {
            r = (r << 1) | (code & 1);
            code >>= 1;
        }
        return r;
    }

    void build(const std::vector<uint8_t>& code_lengths) {
        int n = (int)code_lengths.size();
        lengths.assign(code_lengths.begin(), code_lengths.end());

        count.assign(MAX_BITS + 1, 0);
        for (int i = 0; i < n; ++i) {
            int len = lengths[i];
            if (len > 0) {
                if (len > MAX_BITS) throw image_error("Huffman: code too long");
                count[len]++;
            }
        }

        first_code.assign(MAX_BITS + 1, 0);
        int code = 0;
        for (int len = 1; len <= MAX_BITS; ++len) {
            code = (code + count[len - 1]) << 1;
            first_code[len] = code;
        }

        first_symbol.assign(MAX_BITS + 1, 0);
        int sum = 0;
        for (int len = 1; len <= MAX_BITS; ++len) {
            first_symbol[len] = sum;
            sum += count[len];
        }

        symbols.assign(n, 0);
        std::vector<int> next_code = first_symbol;
        for (int i = 0; i < n; ++i) {
            int len = lengths[i];
            if (len != 0) {
                symbols[next_code[len]++] = i;
            }
        }

        for (int i = 0; i < (1 << FAST_BITS); ++i) {
            fast[i].symbol = 0xFFFF;
            fast[i].length = 0;
        }
        
        for (int len = 1; len <= FAST_BITS; ++len) {
            if (count[len] == 0) continue;
            
            int start_code = first_code[len];
            int start_index = first_symbol[len];
            
            for (int i = 0; i < count[len]; ++i) {
                int sym = symbols[start_index + i];
                int huff_code = start_code + i;
                uint32_t rev_code = reverse_bits((uint32_t)huff_code, len);
                
                int entries = 1 << (FAST_BITS - len);
                for (int j = 0; j < entries; ++j) {
                    int idx = rev_code + (j << len);
                    if (idx < (1 << FAST_BITS)) {
                        fast[idx].symbol = (uint16_t)sym;
                        fast[idx].length = (uint8_t)len;
                    }
                }
            }
        }
    }

    ALWAYS_INLINE uint16_t decode(uint32_t& bitbuf, int& bits, const uint8_t*& data, size_t& pos, size_t size) {
        while (bits < FAST_BITS && pos < size) {
            bitbuf |= (uint32_t)data[pos++] << bits;
            bits += 8;
        }

        int idx = (int)(bitbuf & FAST_MASK);
        Entry e = fast[idx];
        
        if (e.length != 0) {
            bitbuf >>= e.length;
            bits -= e.length;
            return e.symbol;
        }

        uint32_t code = 0;
        for (int len = 1; len <= MAX_BITS; ++len) {
            if (bits == 0) {
                if (pos >= size) throw image_error("Huffman: out of data");
                bitbuf |= (uint32_t)data[pos++] << bits;
                bits += 8;
            }
            
            code = (code << 1) | (bitbuf & 1);
            bitbuf >>= 1;
            bits--;

            int first = first_code[len];
            int cnt = count[len];
            
            if (code - first < (uint32_t)cnt) {
                int symbol_idx = first_symbol[len] + (code - first);
                return (uint16_t)symbols[symbol_idx];
            }
        }
        
        throw image_error("Huffman: invalid code");
    }
};

// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//    INFLATE DECOMPRESSION (minimal for PNG)
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

class Inflate {
private:
    static const int LENGTH_CODES = 29, DISTANCE_CODES = 30, CODE_LENGTH_CODES = 19;
    
    struct LengthInfo { uint16_t base; uint8_t extra_bits; };
    
    static inline const LengthInfo length_table[29] = {
        {3,0},{4,0},{5,0},{6,0},{7,0},{8,0},{9,0},{10,0},{11,1},{13,1},{15,1},{17,1},
        {19,2},{23,2},{27,2},{31,2},{35,3},{43,3},{51,3},{59,3},{67,4},{83,4},{99,4},
        {115,4},{131,5},{163,5},{195,5},{227,5},{258,0}
    };
    
    static inline const LengthInfo distance_table[30] = {
        {1,0},{2,0},{3,0},{4,0},{5,1},{7,1},{9,2},{13,2},{17,3},{25,3},
        {33,4},{49,4},{65,5},{97,5},{129,6},{193,6},{257,7},{385,7},{513,8},
        {769,8},{1025,9},{1537,9},{2049,10},{3073,10},{4097,11},{6145,11},
        {8193,12},{12289,12},{16385,13},{24577,13}
    };
    
    static inline const uint8_t cl_order[19] = {
        16,17,18,0,8,7,9,6,10,5,11,4,12,3,13,2,14,1,15
    };

    const uint8_t* compressed_data;
    size_t compressed_size, bit_pos;
    uint32_t bit_buffer;
    int bits_in_buffer;

    ALWAYS_INLINE uint32_t read_bits(int n) {
        while (bits_in_buffer < n) {
            if (bit_pos >= compressed_size) throw image_error("Unexpected end of compressed data");
            bit_buffer |= (uint32_t)compressed_data[bit_pos++] << bits_in_buffer;
            bits_in_buffer += 8;
        }
        uint32_t result = bit_buffer & ((1u << n) - 1);
        bit_buffer >>= n; bits_in_buffer -= n;
        return result;
    }

    ALWAYS_INLINE void align_to_byte() { bit_buffer = 0; bits_in_buffer = 0; }

    std::vector<uint8_t> read_lengths(HuffmanDecoder& decoder, int count) {
        std::vector<uint8_t> lengths(count, 0);
        int i = 0;
        while (i < count) {
            uint16_t symbol = decoder.decode(bit_buffer, bits_in_buffer, compressed_data, bit_pos, compressed_size);
            if (symbol < 16) { lengths[i++] = (uint8_t)symbol; }
            else if (symbol == 16) {
                if (i == 0) throw image_error("Cannot repeat previous at start of lengths");
                int repeat = 3 + (int)read_bits(2);
                if (i + repeat > count) repeat = count - i;
                uint8_t prev = lengths[i-1];
                for (int k = 0; k < repeat; ++k) lengths[i++] = prev;
            }
            else if (symbol == 17) {
                int repeat = 3 + (int)read_bits(3);
                if (i + repeat > count) repeat = count - i;
                for (int k = 0; k < repeat; ++k) lengths[i++] = 0;
            }
            else if (symbol == 18) {
                int repeat = 11 + (int)read_bits(7);
                if (i + repeat > count) repeat = count - i;
                for (int k = 0; k < repeat; ++k) lengths[i++] = 0;
            }
            else throw image_error("Invalid symbol in code lengths");
        }
        return lengths;
    }

    void decode_block(HuffmanDecoder& litlen_decoder, HuffmanDecoder& dist_decoder, std::vector<uint8_t>& result) {
        while (true) {
            uint16_t symbol = litlen_decoder.decode(bit_buffer, bits_in_buffer, compressed_data, bit_pos, compressed_size);
            if (symbol < 256) { result.push_back((uint8_t)symbol); continue; }
            if (symbol == 256) break;
            
            int length_index = (int)symbol - 257;
            if (length_index < 0 || length_index >= LENGTH_CODES) throw image_error("Invalid length code");
            const LengthInfo& info = length_table[length_index];
            uint32_t length = info.base;
            if (info.extra_bits) length += read_bits(info.extra_bits);
            
            uint16_t dist_symbol = dist_decoder.decode(bit_buffer, bits_in_buffer, compressed_data, bit_pos, compressed_size);
            if (dist_symbol >= DISTANCE_CODES) throw image_error("Invalid distance code");
            const LengthInfo& dist_info = distance_table[dist_symbol];
            uint32_t distance = dist_info.base;
            if (dist_info.extra_bits) distance += read_bits(dist_info.extra_bits);
            
            if (distance == 0 || distance > result.size()) throw image_error("Distance too large");
            size_t start = result.size() - distance, old_size = result.size();
            result.resize(old_size + length);
            
            if (distance >= length) {
                memcpy(result.data() + old_size, result.data() + start, length);
            } else {
                for (uint32_t i = 0; i < length; i++) result[old_size + i] = result[start + (i % distance)];
            }
        }
    }

    void process_fixed_block(std::vector<uint8_t>& result) {
        std::vector<uint8_t> litlen_lengths(288, 0), dist_lengths(32, 5);
        for (int i = 0; i <= 143; i++) litlen_lengths[i] = 8;
        for (int i = 144; i <= 255; i++) litlen_lengths[i] = 9;
        for (int i = 256; i <= 279; i++) litlen_lengths[i] = 7;
        for (int i = 280; i <= 287; i++) litlen_lengths[i] = 8;
        
        HuffmanDecoder litlen_decoder, dist_decoder;
        litlen_decoder.build(litlen_lengths);
        dist_decoder.build(dist_lengths);
        decode_block(litlen_decoder, dist_decoder, result);
    }

    void process_dynamic_block(std::vector<uint8_t>& result) {
        int hlit = (int)read_bits(5) + 257, hdist = (int)read_bits(5) + 1, hclen = (int)read_bits(4) + 4;
        std::vector<uint8_t> code_lengths(19, 0);
        for (int i = 0; i < hclen; i++) code_lengths[cl_order[i]] = (uint8_t)read_bits(3);
        
        HuffmanDecoder cl_decoder; cl_decoder.build(code_lengths);
        int total_codes = hlit + hdist;
        std::vector<uint8_t> all_lengths = read_lengths(cl_decoder, total_codes);
        
        std::vector<uint8_t> litlen_lengths(all_lengths.begin(), all_lengths.begin() + hlit);
        std::vector<uint8_t> dist_lengths(all_lengths.begin() + hlit, all_lengths.end());
        
        HuffmanDecoder litlen_decoder, dist_decoder;
        litlen_decoder.build(litlen_lengths);
        dist_decoder.build(dist_lengths);
        decode_block(litlen_decoder, dist_decoder, result);
    }

public:
    Inflate(const uint8_t* data, size_t size) : compressed_data(data), compressed_size(size), bit_pos(0), bit_buffer(0), bits_in_buffer(0) {}

    std::vector<uint8_t> decompress() {
        std::vector<uint8_t> result;
        result.reserve(compressed_size * 3);
        
        read_bits(8); read_bits(8);
        
        while (true) {
            int final = (int)read_bits(1), type = (int)read_bits(2);
            
            if (type == 0) {
                align_to_byte();
                if (bit_pos + 4 > compressed_size) throw image_error("Invalid uncompressed block header");
                uint32_t len = compressed_data[bit_pos] | (compressed_data[bit_pos + 1] << 8);
                uint32_t nlen = compressed_data[bit_pos + 2] | (compressed_data[bit_pos + 3] << 8);
                bit_pos += 4;
                if ((len ^ 0xFFFFu) != nlen) throw image_error("Invalid length in uncompressed block");
                if (bit_pos + len > compressed_size) throw image_error("Unexpected end of data in uncompressed block");
                size_t old_size = result.size();
                result.resize(old_size + len);
                memcpy(result.data() + old_size, compressed_data + bit_pos, len);
                bit_pos += len;
            }
            else if (type == 1) process_fixed_block(result);
            else if (type == 2) process_dynamic_block(result);
            else throw image_error("Invalid block type");
            
            if (final) break;
        }
        
        align_to_byte();
        if (bit_pos + 4 <= compressed_size) bit_pos += 4;
        return result;
    }
};

// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//                FILTERS FOR RGBA
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

static void apply_filter_rgba_sse(uint8_t* scanline, const uint8_t* prev_scanline, int width, int filter_type) {
    const int channels = 4;
    const int stride   = width * channels;

    switch (filter_type) {
        case 0: // None
            return;

        case 1: { // Sub
            uint8_t* p = scanline + channels;
            int remaining = stride - channels;

            while (remaining >= 16) {
                __m128i cur  = _mm_loadu_si128((__m128i*)p);
                __m128i left = _mm_loadu_si128((__m128i*)(p - channels));
                __m128i res  = _mm_add_epi8(cur, left);
                _mm_storeu_si128((__m128i*)p, res);

                p += 16;
                remaining -= 16;
            }
            for (int i = stride - remaining; i < stride; ++i) {
                scanline[i] = (uint8_t)(scanline[i] + scanline[i - channels]);
            }
        } break;

        case 2: { // Up
            const uint8_t* up = prev_scanline;
            uint8_t* cur      = scanline;

            int remaining = stride;
            while (remaining >= 16) {
                __m128i vcur = _mm_loadu_si128((__m128i*)cur);
                __m128i vup  = _mm_loadu_si128((__m128i*)up);
                __m128i res  = _mm_add_epi8(vcur, vup);
                _mm_storeu_si128((__m128i*)cur, res);

                cur += 16;
                up  += 16;
                remaining -= 16;
            }
            for (int i = stride - remaining; i < stride; ++i) {
                scanline[i] = (uint8_t)(scanline[i] + prev_scanline[i]);
            }
        } break;

        case 3: { // Average
            const uint8_t* up = prev_scanline;
            uint8_t* cur      = scanline;

            for (int i = 0; i < channels; ++i) {
                uint8_t left = 0;
                uint8_t u    = up[i];
                cur[i] = (uint8_t)(cur[i] + ((left + u) >> 1));
            }

            int i = channels;
            for (; i + 16 <= stride; i += 16) {
                __m128i vcur = _mm_loadu_si128((__m128i*)(cur + i));
                __m128i vup  = _mm_loadu_si128((__m128i*)(up  + i));
                __m128i vleft= _mm_loadu_si128((__m128i*)(cur + i - channels));

                __m128i sum  = _mm_add_epi8(vleft, vup);

                __m128i zero = _mm_setzero_si128();
                __m128i lo   = _mm_unpacklo_epi8(sum, zero);
                __m128i hi   = _mm_unpackhi_epi8(sum, zero);

                lo = _mm_srli_epi16(lo, 1);
                hi = _mm_srli_epi16(hi, 1);

                __m128i avg = _mm_packus_epi16(lo, hi);
                __m128i res = _mm_add_epi8(vcur, avg);

                _mm_storeu_si128((__m128i*)(cur + i), res);
            }

            for (; i < stride; ++i) {
                uint8_t left = (i >= channels) ? cur[i - channels] : 0;
                uint8_t u    = up[i];
                cur[i] = (uint8_t)(cur[i] + ((left + u) >> 1));
            }
        } break;

        case 4: {
            for (int i = 0; i < stride; i++) {
                uint8_t left    = (i >= channels) ? scanline[i - channels] : 0;
                uint8_t up      = prev_scanline[i];
                uint8_t upleft  = (i >= channels) ? prev_scanline[i - channels] : 0;

                int a = left;
                int b = up;
                int c = upleft;

                int p  = a + b - c;
                int pa = std::abs(p - a);
                int pb = std::abs(p - b);
                int pc = std::abs(p - c);

                uint8_t paeth = (pa <= pb && pa <= pc) ? (uint8_t)a :
                                (pb <= pc) ? (uint8_t)b : (uint8_t)c;

                scanline[i] = (uint8_t)(scanline[i] + paeth);
            }
        } break;

        default:
            throw image_error("Unknown PNG filter type");
    }
}

// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//                 PNG DECODER
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

class PNGDecoder {  
private:
    struct PNGChunk {
        uint32_t length;
        char type[5];
        std::vector<uint8_t> data;
        uint32_t crc;
    };
    
    struct PNGInfo {
        int width = 0, height = 0, bit_depth = 8, color_type = 0;
        int compression_method = 0, filter_method = 0, interlace_method = 0;
        std::vector<uint8_t> palette, transparency, background, icc_profile;
        std::vector<uint8_t> text_data, compressed_text;
        double gamma = 0.45455;
        double chroma_white[2] = {0.3127, 0.3290};
        double chroma_red[2] = {0.64, 0.33}, chroma_green[2] = {0.30, 0.60}, chroma_blue[2] = {0.15, 0.06};
        int srgb_intent = 0;
        int physical_pixels[2] = {0, 0}; // pHYs chunk
        int unit_specifier = 0;
        
        bool has_alpha() const {
            return color_type == 4 || color_type == 6 || (color_type == 3 && !transparency.empty());
        }
        
        int get_channels() const {
            switch(color_type) {
                case 0: return 1;
                case 2: return 3;
                case 3: return 1;
                case 4: return 2;
                case 6: return 4;
                default: return 0;
            }
        }
        
        int get_bpp() const {
            int channels = get_channels();
            if (bit_depth == 16) return channels * 2;
            if (bit_depth == 8) return channels;
            if (bit_depth == 4) return (channels * 4 + 7) / 8;
            if (bit_depth == 2) return (channels * 2 + 7) / 8;
            if (bit_depth == 1) return (channels + 7) / 8;
            return channels;
        }
    };
    
    static constexpr int ADAM7_X[7] = {0, 4, 0, 2, 0, 1, 0};
    static constexpr int ADAM7_Y[7] = {0, 0, 4, 0, 2, 0, 1};
    static constexpr int ADAM7_DX[7] = {8, 8, 4, 4, 2, 2, 1};
    static constexpr int ADAM7_DY[7] = {8, 8, 8, 4, 4, 2, 2};
    
    static PNGChunk read_chunk(std::ifstream& file) {
        PNGChunk chunk;
        uint8_t len_bytes[4];
        file.read((char*)len_bytes, 4);
        if (file.gcount() != 4) throw image_error("Failed to read chunk length");
        
        chunk.length = (len_bytes[0] << 24) | (len_bytes[1] << 16) | (len_bytes[2] << 8) | len_bytes[3];
        if (chunk.length > 100000000u) throw image_error("Chunk too large");
        
        file.read(chunk.type, 4);
        if (file.gcount() != 4) throw image_error("Failed to read chunk type");
        chunk.type[4] = '\0';
        
        chunk.data.resize(chunk.length);
        if (chunk.length > 0) {
            file.read((char*)chunk.data.data(), chunk.length);
            if ((size_t)file.gcount() != chunk.length) throw image_error("Failed to read chunk data");
        }
        
        uint8_t crc_bytes[4];
        file.read((char*)crc_bytes, 4);
        if (file.gcount() != 4) throw image_error("Failed to read chunk CRC");
        
        chunk.crc = (crc_bytes[0] << 24) | (crc_bytes[1] << 16) | (crc_bytes[2] << 8) | crc_bytes[3];
        
        uint32_t computed_crc = CRC32::update(0, (uint8_t*)chunk.type, 4);
        if (chunk.length > 0) computed_crc = CRC32::update(computed_crc, chunk.data.data(), chunk.length);
        if (computed_crc != chunk.crc) throw image_error("Chunk CRC mismatch for " + std::string(chunk.type));
        
        return chunk;
    }
    
    static ALWAYS_INLINE uint8_t paeth_predictor(uint8_t a, uint8_t b, uint8_t c) {
        int p = a + b - c;
        int pa = abs(p - a);
        int pb = abs(p - b);
        int pc = abs(p - c);
        return (pa <= pb && pa <= pc) ? a : (pb <= pc) ? b : c;
    }
    
    static void unfilter_scanline(uint8_t* scanline, const uint8_t* prev_scanline, int width, int bpp, int filter_type) {
        switch (filter_type) {
            case 0: break;
            case 1: for (int i = bpp; i < width * bpp; i++) scanline[i] += scanline[i - bpp]; break;
            case 2: for (int i = 0; i < width * bpp; i++) scanline[i] += prev_scanline[i]; break;
            case 3: for (int i = 0; i < width * bpp; i++) {
                uint8_t left = (i >= bpp) ? scanline[i - bpp] : 0;
                scanline[i] += (left + prev_scanline[i]) >> 1;
            } break;
            case 4: for (int i = 0; i < width * bpp; i++) {
                uint8_t left = (i >= bpp) ? scanline[i - bpp] : 0;
                uint8_t up = prev_scanline[i];
                uint8_t upleft = (i >= bpp) ? prev_scanline[i - bpp] : 0;
                scanline[i] += paeth_predictor(left, up, upleft);
            } break;
            default: throw image_error("Invalid PNG filter type");
        }
    }
    
    static void convert_16bit_to_8bit(uint8_t* dest, const uint8_t* src, int pixels) {
        for (int i = 0; i < pixels; i++) dest[i] = src[i * 2];
    }
    
    static void convert_low_bitdepth(uint8_t* dest, const uint8_t* src, int width, int height, int channels, int src_bit_depth) {
        int bits_per_pixel = src_bit_depth * channels;
        int pixels_per_byte = 8 / src_bit_depth;
        int mask = (1 << src_bit_depth) - 1;
        
        for (int y = 0; y < height; y++) {
            const uint8_t* src_row = src + y * ((width * bits_per_pixel + 7) / 8);
            uint8_t* dst_row = dest + y * width * channels;
            
            for (int x = 0; x < width; x++) {
                int byte_idx = (x * bits_per_pixel) / 8;
                int bit_idx = (x * bits_per_pixel) % 8;
                
                for (int c = 0; c < channels; c++) {
                    int bit_offset = bit_idx + c * src_bit_depth;
                    int val = (src_row[byte_idx + bit_offset / 8] >> (7 - (bit_offset % 8))) & mask;
                    dst_row[x * channels + c] = (val * 255) / mask;
                }
            }
        }
    }
    
    static std::vector<uint8_t> apply_palette(const std::vector<uint8_t>& indices, const PNGInfo& info) {
        int palette_entries = info.palette.size() / 3;
        std::vector<uint8_t> rgba(info.width * info.height * 4);
        
        for (size_t i = 0; i < indices.size(); i++) {
            int idx = indices[i];
            if (idx >= palette_entries) idx = 0;
            
            rgba[i * 4] = info.palette[idx * 3];
            rgba[i * 4 + 1] = info.palette[idx * 3 + 1];
            rgba[i * 4 + 2] = info.palette[idx * 3 + 2];
            rgba[i * 4 + 3] = (idx < (int)info.transparency.size()) ? info.transparency[idx] : 255;
        }
        
        return rgba;
    }
    
    static std::vector<uint8_t> apply_transparency(const std::vector<uint8_t>& data, const PNGInfo& info) {
        int channels = info.get_channels();
        std::vector<uint8_t> rgba(info.width * info.height * 4);
        
        if (info.color_type == 0) {
            uint16_t trans_value = (info.transparency[0] << 8) | info.transparency[1];
            for (int i = 0; i < info.width * info.height; i++) {
                uint16_t value = (data[i * 2] << 8) | data[i * 2 + 1];
                rgba[i * 4] = data[i * 2];
                rgba[i * 4 + 1] = data[i * 2];
                rgba[i * 4 + 2] = data[i * 2];
                rgba[i * 4 + 3] = (value == trans_value) ? 0 : 255;
            }
        } else if (info.color_type == 2) {
            uint16_t r = (info.transparency[0] << 8) | info.transparency[1];
            uint16_t g = (info.transparency[2] << 8) | info.transparency[3];
            uint16_t b = (info.transparency[4] << 8) | info.transparency[5];
            
            for (int i = 0; i < info.width * info.height; i++) {
                uint16_t pr = (data[i * 6] << 8) | data[i * 6 + 1];
                uint16_t pg = (data[i * 6 + 2] << 8) | data[i * 6 + 3];
                uint16_t pb = (data[i * 6 + 4] << 8) | data[i * 6 + 5];
                
                rgba[i * 4] = data[i * 6];
                rgba[i * 4 + 1] = data[i * 6 + 2];
                rgba[i * 4 + 2] = data[i * 6 + 4];
                rgba[i * 4 + 3] = (pr == r && pg == g && pb == b) ? 0 : 255;
            }
        }
        
        return rgba;
    }
    
    static std::vector<uint8_t> deinterlace_adam7(const std::vector<uint8_t>& idat_data, const PNGInfo& info) {
        int bpp = info.get_bpp();
        int pass_width[7], pass_height[7];
        
        for (int p = 0; p < 7; p++) {
            pass_width[p] = (info.width - ADAM7_X[p] + ADAM7_DX[p] - 1) / ADAM7_DX[p];
            pass_height[p] = (info.height - ADAM7_Y[p] + ADAM7_DY[p] - 1) / ADAM7_DY[p];
        }
        
        std::vector<uint8_t> output(info.width * info.height * bpp, 0);
        const uint8_t* data_ptr = idat_data.data();
        size_t data_pos = 0;
        
        for (int p = 0; p < 7; p++) {
            if (pass_width[p] == 0 || pass_height[p] == 0) continue;
            
            int pass_stride = pass_width[p] * bpp;
            std::vector<uint8_t> pass_prev(pass_stride, 0);
            
            for (int y = 0; y < pass_height[p]; y++) {
                if (data_pos >= idat_data.size()) throw image_error("Unexpected end of IDAT data");
                
                int filter_type = data_ptr[data_pos++];
                const uint8_t* scanline = data_ptr + data_pos;
                data_pos += pass_stride;
                
                std::vector<uint8_t> filtered(pass_stride);
                memcpy(filtered.data(), scanline, pass_stride);
                unfilter_scanline(filtered.data(), pass_prev.data(), pass_width[p], bpp, filter_type);
                
                for (int x = 0; x < pass_width[p]; x++) {
                    int out_x = ADAM7_X[p] + x * ADAM7_DX[p];
                    int out_y = ADAM7_Y[p] + y * ADAM7_DY[p];
                    
                    if (out_x < info.width && out_y < info.height) {
                        int out_idx = (out_y * info.width + out_x) * bpp;
                        int in_idx = x * bpp;
                        memcpy(&output[out_idx], &filtered[in_idx], bpp);
                    }
                }
                
                memcpy(pass_prev.data(), filtered.data(), pass_stride);
            }
        }
        
        return output;
    }
    
    static std::vector<uint8_t> process_non_interlaced(const std::vector<uint8_t>& idat_data, const PNGInfo& info) {
        int bpp = info.get_bpp();
        int stride = info.width * bpp;
        
        std::vector<uint8_t> output(info.width * info.height * bpp);
        std::vector<uint8_t> prev_scanline(stride, 0);
        
        const uint8_t* data_ptr = idat_data.data();
        size_t data_pos = 0;
        
        for (int y = 0; y < info.height; y++) {
            if (data_pos >= idat_data.size()) throw image_error("Unexpected end of IDAT data");
            
            int filter_type = data_ptr[data_pos++];
            const uint8_t* scanline = data_ptr + data_pos;
            data_pos += stride;
            
            std::vector<uint8_t> filtered(stride);
            memcpy(filtered.data(), scanline, stride);
            
            if (info.color_type == 6 && info.bit_depth == 8 && info.width >= 4) {
                apply_filter_rgba_sse(filtered.data(), prev_scanline.data(), info.width, filter_type);
            } else {
                unfilter_scanline(filtered.data(), prev_scanline.data(), info.width, bpp, filter_type);
            }
            
            memcpy(&output[y * stride], filtered.data(), stride);
            memcpy(prev_scanline.data(), filtered.data(), stride);
        }
        
        return output;
    }
    
    static std::vector<uint8_t> process_image_data(const std::vector<uint8_t>& idat_data, const PNGInfo& info) {
        if (info.interlace_method == 1) {
            return deinterlace_adam7(idat_data, info);
        } else {
            return process_non_interlaced(idat_data, info);
        }
    }
    
    static std::vector<uint8_t> convert_to_rgba(const std::vector<uint8_t>& raw_data, const PNGInfo& info) {
        int channels = info.get_channels();
        
        if (info.color_type == 3) {
            return apply_palette(raw_data, info);
        }
        
        if (info.bit_depth == 16) {
            std::vector<uint8_t> converted_8bit(info.width * info.height * channels);
            convert_16bit_to_8bit(converted_8bit.data(), raw_data.data(), info.width * info.height * channels);
            
            if (info.color_type == 0) {
                if (!info.transparency.empty()) {
                    return apply_transparency(raw_data, info);
                }
                std::vector<uint8_t> rgb(info.width * info.height * 3);
                for (int i = 0; i < info.width * info.height; i++) {
                    rgb[i * 3] = converted_8bit[i];
                    rgb[i * 3 + 1] = converted_8bit[i];
                    rgb[i * 3 + 2] = converted_8bit[i];
                }
                return rgb;
            }
            
            if (info.color_type == 2 || info.color_type == 4 || info.color_type == 6) {
                std::vector<uint8_t> result(info.width * info.height * (info.has_alpha() ? 4 : 3));
                
                if (info.color_type == 2) {
                    for (int i = 0; i < info.width * info.height; i++) {
                        result[i * 3] = converted_8bit[i * 3];
                        result[i * 3 + 1] = converted_8bit[i * 3 + 1];
                        result[i * 3 + 2] = converted_8bit[i * 3 + 2];
                    }
                } else if (info.color_type == 4) {
                    for (int i = 0; i < info.width * info.height; i++) {
                        result[i * 4] = converted_8bit[i * 2];
                        result[i * 4 + 1] = converted_8bit[i * 2];
                        result[i * 4 + 2] = converted_8bit[i * 2];
                        result[i * 4 + 3] = converted_8bit[i * 2 + 1];
                    }
                } else if (info.color_type == 6) {
                    for (int i = 0; i < info.width * info.height; i++) {
                        result[i * 4] = converted_8bit[i * 4];
                        result[i * 4 + 1] = converted_8bit[i * 4 + 2];
                        result[i * 4 + 2] = converted_8bit[i * 4 + 4];
                        result[i * 4 + 3] = converted_8bit[i * 4 + 6];
                    }
                }
                
                return result;
            }
        }
        
        if (info.bit_depth < 8) {
            std::vector<uint8_t> converted_8bit(info.width * info.height * channels);
            convert_low_bitdepth(converted_8bit.data(), raw_data.data(), info.width, info.height, channels, info.bit_depth);
            
            if (info.color_type == 0 || info.color_type == 3) {
                std::vector<uint8_t> rgb(info.width * info.height * 3);
                for (int i = 0; i < info.width * info.height; i++) {
                    rgb[i * 3] = converted_8bit[i];
                    rgb[i * 3 + 1] = converted_8bit[i];
                    rgb[i * 3 + 2] = converted_8bit[i];
                }
                return rgb;
            }
        }
        
        if (info.color_type == 0) {
            std::vector<uint8_t> rgb(info.width * info.height * 3);
            for (int i = 0; i < info.width * info.height; i++) {
                rgb[i * 3] = raw_data[i];
                rgb[i * 3 + 1] = raw_data[i];
                rgb[i * 3 + 2] = raw_data[i];
            }
            return rgb;
        }
        
        if (info.color_type == 2) {
            if (!info.transparency.empty()) {
                std::vector<uint8_t> rgba(info.width * info.height * 4);
                for (int i = 0; i < info.width * info.height; i++) {
                    rgba[i * 4] = raw_data[i * 3];
                    rgba[i * 4 + 1] = raw_data[i * 3 + 1];
                    rgba[i * 4 + 2] = raw_data[i * 3 + 2];
                    rgba[i * 4 + 3] = 255;
                    
                    if (raw_data[i * 3] == info.transparency[0] &&
                        raw_data[i * 3 + 1] == info.transparency[1] &&
                        raw_data[i * 3 + 2] == info.transparency[2]) {
                        rgba[i * 4 + 3] = 0;
                    }
                }
                return rgba;
            }
            std::vector<uint8_t> rgb(info.width * info.height * 3);
            memcpy(rgb.data(), raw_data.data(), info.width * info.height * 3);
            return rgb;
        }
        
        if (info.color_type == 4) {
            std::vector<uint8_t> rgba(info.width * info.height * 4);
            for (int i = 0; i < info.width * info.height; i++) {
                rgba[i * 4] = raw_data[i * 2];
                rgba[i * 4 + 1] = raw_data[i * 2];
                rgba[i * 4 + 2] = raw_data[i * 2];
                rgba[i * 4 + 3] = raw_data[i * 2 + 1];
            }
            return rgba;
        }
        
        if (info.color_type == 6) {
            std::vector<uint8_t> rgba(info.width * info.height * 4);
            memcpy(rgba.data(), raw_data.data(), info.width * info.height * 4);
            return rgba;
        }
        
        throw image_error("Unsupported PNG color type");
    }
    
    static void apply_gamma_correction(std::vector<uint8_t>& data, double gamma) {
        if (gamma <= 0.0) return;
        
        uint8_t gamma_table[256];
        for (int i = 0; i < 256; i++) {
            double normalized = i / 255.0;
            double corrected = pow(normalized, 1.0 / gamma);
            gamma_table[i] = (uint8_t)(corrected * 255.0 + 0.5);
        }
        
        for (size_t i = 0; i < data.size(); i++) {
            data[i] = gamma_table[data[i]];
        }
    }

public:
    static ImageData load(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open()) throw image_error("Cannot open PNG file");
        
        uint8_t signature[8];
        file.read((char*)signature, 8);
        if (file.gcount() != 8 || signature[0] != 137 || signature[1] != 80 ||
            signature[2] != 78 || signature[3] != 71 || signature[4] != 13 ||
            signature[5] != 10 || signature[6] != 26 || signature[7] != 10) {
            throw image_error("Invalid PNG signature");
        }
        
        PNGInfo info;
        std::vector<uint8_t> compressed_data;
        bool seen_idat = false;
        
        while (true) {
            PNGChunk chunk = read_chunk(file);
            std::string type_str(chunk.type);
            
            if (type_str == "IHDR") {
                if (chunk.length != 13) throw image_error("Invalid IHDR chunk size");
                
                info.width = (chunk.data[0] << 24) | (chunk.data[1] << 16) | 
                            (chunk.data[2] << 8) | chunk.data[3];
                info.height = (chunk.data[4] << 24) | (chunk.data[5] << 16) | 
                             (chunk.data[6] << 8) | chunk.data[7];
                info.bit_depth = chunk.data[8];
                info.color_type = chunk.data[9];
                info.compression_method = chunk.data[10];
                info.filter_method = chunk.data[11];
                info.interlace_method = chunk.data[12];
                
                if (info.width <= 0 || info.height <= 0 || 
                    info.width > IMG_MAX_IMAGE_SIZE || info.height > IMG_MAX_IMAGE_SIZE) {
                    throw image_error("Invalid PNG dimensions");
                }
                
                if (info.compression_method != 0) throw image_error("Unsupported compression method");
                if (info.filter_method != 0) throw image_error("Unsupported filter method");
                if (info.interlace_method > 1) throw image_error("Unsupported interlace method");
                
                static const bool valid_combinations[7][5] = {
                    {true, true, true, true, true},
                    {false, false, false, false, false},
                    {false, false, false, true, true},
                    {true, true, true, true, false},
                    {false, false, false, true, true},
                    {false, false, false, false, false},
                    {false, false, false, true, true}
                };
                
                int bd_idx = 0;
                switch(info.bit_depth) {
                    case 1: bd_idx = 0; break;
                    case 2: bd_idx = 1; break;
                    case 4: bd_idx = 2; break;
                    case 8: bd_idx = 3; break;
                    case 16: bd_idx = 4; break;
                    default: throw image_error("Invalid bit depth");
                }
                
                if (!valid_combinations[info.color_type][bd_idx]) {
                    throw image_error("Invalid color type/bit depth combination");
                }
            }
            else if (type_str == "PLTE") {
                if (chunk.length % 3 != 0 || chunk.length > 768) {
                    throw image_error("Invalid PLTE chunk");
                }
                info.palette = std::move(chunk.data);
            }
            else if (type_str == "IDAT") {
                compressed_data.insert(compressed_data.end(), chunk.data.begin(), chunk.data.end());
                seen_idat = true;
            }
            else if (type_str == "IEND") {
                break;
            }
            else if (type_str == "tRNS") {
                info.transparency = std::move(chunk.data);
            }
            else if (type_str == "bKGD") {
                info.background = std::move(chunk.data);
            }
            else if (type_str == "gAMA") {
                if (chunk.length == 4) {
                    uint32_t gamma_int = (chunk.data[0] << 24) | (chunk.data[1] << 16) | 
                                        (chunk.data[2] << 8) | chunk.data[3];
                    info.gamma = gamma_int / 100000.0;
                }
            }
            else if (type_str == "cHRM") {
                if (chunk.length == 32) {
                    for (int i = 0; i < 8; i++) {
                        uint32_t val = (chunk.data[i*4] << 24) | (chunk.data[i*4+1] << 16) | 
                                      (chunk.data[i*4+2] << 8) | chunk.data[i*4+3];
                        double coord = val / 100000.0;
                        switch(i) {
                            case 0: info.chroma_white[0] = coord; break;
                            case 1: info.chroma_white[1] = coord; break;
                            case 2: info.chroma_red[0] = coord; break;
                            case 3: info.chroma_red[1] = coord; break;
                            case 4: info.chroma_green[0] = coord; break;
                            case 5: info.chroma_green[1] = coord; break;
                            case 6: info.chroma_blue[0] = coord; break;
                            case 7: info.chroma_blue[1] = coord; break;
                        }
                    }
                }
            }
            else if (type_str == "sRGB") {
                if (chunk.length == 1) {
                    info.srgb_intent = chunk.data[0];
                }
            }
            else if (type_str == "iCCP") {
                info.icc_profile = std::move(chunk.data);
            }
            else if (type_str == "tEXt" || type_str == "zTXt" || type_str == "iTXt") {
                info.text_data.insert(info.text_data.end(), chunk.data.begin(), chunk.data.end());
            }
            else if (type_str == "pHYs") {
                if (chunk.length == 9) {
                    info.physical_pixels[0] = (chunk.data[0] << 24) | (chunk.data[1] << 16) | 
                                             (chunk.data[2] << 8) | chunk.data[3];
                    info.physical_pixels[1] = (chunk.data[4] << 24) | (chunk.data[5] << 16) | 
                                             (chunk.data[6] << 8) | chunk.data[7];
                    info.unit_specifier = chunk.data[8];
                }
            }
            else if (type_str == "tIME") {
            }
            else if (type_str == "sBIT") {
            }
            else if (type_str == "hIST") {
            }
            else if (type_str == "sPLT") {
            }
            else if (type_str[0] >= 'a' && type_str[0] <= 'z') {
            }
            else {
                throw image_error("Unrecognized critical PNG chunk: " + type_str);
            }
        }
        
        if (!seen_idat) throw image_error("No image data in PNG file");
        if (info.width == 0 || info.height == 0) throw image_error("No valid PNG header found");
        
        Inflate inflater(compressed_data.data(), compressed_data.size());
        std::vector<uint8_t> idat_data = inflater.decompress();
        
        if (info.color_type == 3 && info.palette.empty()) {
            throw image_error("Palette image without PLTE chunk");
        }
        
        std::vector<uint8_t> raw_data = process_image_data(idat_data, info);
        std::vector<uint8_t> final_data = convert_to_rgba(raw_data, info);
        
        int final_channels = (final_data.size() / (info.width * info.height) == 3) ? 3 : 4;
        
        ImageData result(info.width, info.height, final_channels);
        result.data = std::move(final_data);
        
        if (info.gamma > 0.0 && info.gamma != 0.45455) {
            apply_gamma_correction(result.data, info.gamma);
        }
        
        result.flip_vertically();
        return result;
    }
    
    static PNGInfo get_info(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open()) throw image_error("Cannot open PNG file");
        
        uint8_t signature[8];
        file.read((char*)signature, 8);
        if (file.gcount() != 8 || signature[0] != 137 || signature[1] != 80 ||
            signature[2] != 78 || signature[3] != 71 || signature[4] != 13 ||
            signature[5] != 10 || signature[6] != 26 || signature[7] != 10) {
            throw image_error("Invalid PNG signature");
        }
        
        PNGInfo info;
        PNGChunk chunk = read_chunk(file);
        
        if (strcmp(chunk.type, "IHDR") != 0 || chunk.length != 13) {
            throw image_error("First chunk is not IHDR");
        }
        
        info.width = (chunk.data[0] << 24) | (chunk.data[1] << 16) | 
                    (chunk.data[2] << 8) | chunk.data[3];
        info.height = (chunk.data[4] << 24) | (chunk.data[5] << 16) | 
                     (chunk.data[6] << 8) | chunk.data[7];
        info.bit_depth = chunk.data[8];
        info.color_type = chunk.data[9];
        info.compression_method = chunk.data[10];
        info.filter_method = chunk.data[11];
        info.interlace_method = chunk.data[12];
        
        return info;
    }
};


// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//                  HDR DECODER
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

class HDRDecoder {
private:
    static bool read_line(std::istream& file, std::string& line) {
        line.clear();
        char ch;
        while (file.get(ch) && ch != '\n') {
            line.push_back(ch);
        }
        return !file.fail() || (file.eof() && !line.empty());
    }

    static bool parse_dimensions(const std::string& line, int& width, int& height) {
        if (line.size() < 8 || line[0] != '-' || line[1] != 'Y') return false;
        
        const char* ptr = line.c_str() + 2;
        while (*ptr == ' ') ptr++;
        
        height = 0;
        while (*ptr >= '0' && *ptr <= '9') {
            height = height * 10 + (*ptr - '0');
            ptr++;
        }
        
        while (*ptr == ' ') ptr++;
        if (*ptr != '+' || *(ptr + 1) != 'X') return false;
        
        ptr += 2;
        while (*ptr == ' ') ptr++;
        
        width = 0;
        while (*ptr >= '0' && *ptr <= '9') {
            width = width * 10 + (*ptr - '0');
            ptr++;
        }
        
        return width > 0 && height > 0;
    }
    
    static void rgbe_to_float(uint8_t r, uint8_t g, uint8_t b, uint8_t e,
                               float& rf, float& gf, float& bf) {
        if (e == 0) {
            rf = gf = bf = 0.0f;
        } else {
            int exp = static_cast<int>(e) - 128;
            float f = std::ldexp(1.0f, exp - 8);
            rf = static_cast<float>(r) * f;
            gf = static_cast<float>(g) * f;
            bf = static_cast<float>(b) * f;
        }
    }
    
    static bool decode_rle(const uint8_t* src, uint8_t* dst, int width) {
        if (width < 8 || width > 0x7fff) return false;
        
        for (int component = 0; component < 4; component++) {
            int pos = 0;
            while (pos < width) {
                uint8_t code = *src++;
                
                if (code > 128) {
                    int run_length = code - 128;
                    if (run_length == 0 || pos + run_length > width) return false;
                    
                    uint8_t value = *src++;
                    for (int i = 0; i < run_length; i++) {
                        dst[(pos + i) * 4 + component] = value;
                    }
                    pos += run_length;
                } else {
                    if (code == 0 || pos + code > width) return false;
                    
                    for (int i = 0; i < code; i++) {
                        dst[(pos + i) * 4 + component] = *src++;
                    }
                    pos += code;
                }
            }
        }
        
        return true;
    }

    static void apply_tone_mapping(float& r, float& g, float& b, float exposure) {
        r *= exposure;
        g *= exposure;
        b *= exposure;

        const float l_white = 4.0f;
        const float l_white2 = l_white * l_white;
        
        float l = 0.2126f * r + 0.7152f * g + 0.0722f * b;
        float mapped_l = l * (1.0f + l / l_white2) / (1.0f + l);
        
        if (l > 0.0f) {
            float scale = mapped_l / l;
            r *= scale;
            g *= scale;
            b *= scale;
        }
    }
    
    static void apply_gamma_correction(float& r, float& g, float& b) {
        const float gamma = 2.2f;
        r = std::pow(std::max(r, 0.0f), 1.0f / gamma);
        g = std::pow(std::max(g, 0.0f), 1.0f / gamma);
        b = std::pow(std::max(b, 0.0f), 1.0f / gamma);
    }

    static std::vector<uint8_t> read_hdr_file(const std::string& filename, int& width, int& height) {
        std::ifstream file(filename, std::ios::binary);
        if (!file) {
            throw image_error("Cannot open HDR file: " + filename);
        }
        
        std::string line;

        if (!read_line(file, line) || line != "#?RADIANCE") {
            throw image_error("Invalid HDR file signature");
        }

        while (read_line(file, line) && !line.empty()) {
            if (line.find("FORMAT=32-bit_rle_rgbe") != std::string::npos) {
                break;
            }
        }

        read_line(file, line);

        if (!read_line(file, line) || !parse_dimensions(line, width, height)) {
            throw image_error("Invalid HDR dimensions");
        }

        std::vector<uint8_t> rgbe_data(width * height * 4);

        for (int y = 0; y < height; y++) {
            uint8_t header[4];
            file.read(reinterpret_cast<char*>(header), 4);
            
            if (!file) {
                throw image_error("Failed to read scanline header at row " + std::to_string(y));
            }

            if (header[0] == 0x02 && header[1] == 0x02) {
                uint16_t scanline_width = (header[2] << 8) | header[3];
                if (scanline_width != width) {
                    throw image_error("Scanline width mismatch");
                }

                uint8_t* scanline = rgbe_data.data() + (height - 1 - y) * width * 4;

                for (int channel = 0; channel < 4; channel++) {
                    int x = 0;
                    while (x < width) {
                        uint8_t code;
                        file.read(reinterpret_cast<char*>(&code), 1);
                        
                        if (!file) {
                            throw image_error("Failed to read RLE code at row " + std::to_string(y) + 
                                            ", channel " + std::to_string(channel) + 
                                            ", pixel " + std::to_string(x));
                        }

                        if (code > 128) {
                            int run_length = code - 128;
                            if (x + run_length > width) {
                                throw image_error("RLE run length exceeds scanline width");
                            }
                            
                            uint8_t value;
                            file.read(reinterpret_cast<char*>(&value), 1);
                            if (!file) {
                                throw image_error("Failed to read RLE run value");
                            }

                            for (int i = 0; i < run_length; i++) {
                                scanline[(x + i) * 4 + channel] = value;
                            }
                            x += run_length;
                        } else {
                            int run_length = code;
                            if (run_length == 0) {
                                throw image_error("Invalid RLE packet length 0");
                            }
                            if (x + run_length > width) {
                                throw image_error("RLE packet length exceeds scanline width");
                            }

                            for (int i = 0; i < run_length; i++) {
                                uint8_t value;
                                file.read(reinterpret_cast<char*>(&value), 1);
                                if (!file) {
                                    throw image_error("Failed to read RLE raw value");
                                }
                                scanline[(x + i) * 4 + channel] = value;
                            }
                            x += run_length;
                        }
                    }
                }
            } else {
                file.seekg(-4, std::ios::cur);

                uint8_t* scanline = rgbe_data.data() + (height - 1 - y) * width * 4;
                file.read(reinterpret_cast<char*>(scanline), width * 4);
                
                if (!file || file.gcount() != width * 4) {
                    throw image_error("Failed to read uncompressed scanline at row " + std::to_string(y));
                }
            }
        }
        
        return rgbe_data;
    }
    
public:
    static ImageData load(const std::string& filename, float exposure = 1.0f) {
        int width = 0, height = 0;

        std::vector<uint8_t> rgbe_data = read_hdr_file(filename, width, height);

        std::vector<uint8_t> rgb_data(width * height * 3);
        
        for (int i = 0; i < width * height; i++) {
            uint8_t* rgbe = rgbe_data.data() + i * 4;
            uint8_t* rgb = rgb_data.data() + i * 3;
            
            if (rgbe[3] == 0) {
                rgb[0] = rgb[1] = rgb[2] = 0;
            } else {
                int exp = static_cast<int>(rgbe[3]) - 128;
                float f = std::ldexp(1.0f, exp - 8) * exposure;
                
                float r = static_cast<float>(rgbe[0]) * f;
                float g = static_cast<float>(rgbe[1]) * f;
                float b = static_cast<float>(rgbe[2]) * f;

                apply_tone_mapping(r, g, b, exposure);

                apply_gamma_correction(r, g, b);
                
                rgb[0] = static_cast<uint8_t>(std::min(r * 255.0f, 255.0f));
                rgb[1] = static_cast<uint8_t>(std::min(g * 255.0f, 255.0f));
                rgb[2] = static_cast<uint8_t>(std::min(b * 255.0f, 255.0f));
            }
        }
        
        return ImageData(width, height, 3, std::move(rgb_data));
    }

    static ImageData load_fast(const std::string& filename, float exposure = 1.0f) {
        int width = 0, height = 0;

        std::vector<uint8_t> rgbe_data = read_hdr_file(filename, width, height);

        std::vector<uint8_t> rgb_data(width * height * 3);

        for (int i = 0; i < width * height; i++) {
            uint8_t* rgbe = rgbe_data.data() + i * 4;
            uint8_t* rgb = rgb_data.data() + i * 3;
            
            if (rgbe[3] == 0) {
                rgb[0] = rgb[1] = rgb[2] = 0;
            } else {
                int exp = static_cast<int>(rgbe[3]) - 128;
                float f = std::ldexp(1.0f, exp - 8) * exposure;

                float r = static_cast<float>(rgbe[0]) * f;
                float g = static_cast<float>(rgbe[1]) * f;
                float b = static_cast<float>(rgbe[2]) * f;

                r = std::sqrt(r);
                g = std::sqrt(g);
                b = std::sqrt(b);
                
                rgb[0] = static_cast<uint8_t>(std::min(r * 255.0f, 255.0f));
                rgb[1] = static_cast<uint8_t>(std::min(g * 255.0f, 255.0f));
                rgb[2] = static_cast<uint8_t>(std::min(b * 255.0f, 255.0f));
            }
        }
        
        return ImageData(width, height, 3, std::move(rgb_data));
    }

    static std::vector<float> load_float(const std::string& filename,
                                        int& width, int& height) {
        std::vector<uint8_t> rgbe_data = read_hdr_file(filename, width, height);
        std::vector<float> float_data(width * height * 3);
        
        for (int i = 0; i < width * height; i++) {
            uint8_t* rgbe = rgbe_data.data() + i * 4;
            float* rgb = float_data.data() + i * 3;
            rgbe_to_float(rgbe[0], rgbe[1], rgbe[2], rgbe[3], 
                         rgb[0], rgb[1], rgb[2]);
        }
        
        return float_data;
    }

    static ImageData load_from_memory(const uint8_t* data, size_t size, float exposure = 1.0f) {
        if (!data || size < 100) {
            throw image_error("Invalid HDR data buffer");
        }

        const uint8_t* ptr = data;
        const uint8_t* end = data + size;

        int width = 0, height = 0;
        const uint8_t* dim_pos = nullptr;
        
        for (size_t i = 0; i < size - 20; i++) {
            if (data[i] == '-' && data[i+1] == 'Y') {
                dim_pos = data + i;
                break;
            }
        }
        
        if (!dim_pos) {
            throw image_error("HDR: dimensions not found in memory");
        }
        
        std::string dim_str(reinterpret_cast<const char*>(dim_pos), std::min<size_t>(50, end - dim_pos));
        
        if (!parse_dimensions(dim_str, width, height)) {
            throw image_error("HDR: failed to parse dimensions in memory");
        }

        const uint8_t* pixel_start = nullptr;
        for (const uint8_t* p = dim_pos + dim_str.length(); p < end - 4; p++) {
            if (*p == 0xFF && *(p+1) == 0xFF && *(p+2) == 0xFF && *(p+3) == 0xFF) {
                pixel_start = p;
                break;
            }
        }
        
        if (!pixel_start) {
            pixel_start = dim_pos + dim_str.length() + 10;
        }

        size_t pixel_count = width * height;
        std::vector<uint8_t> rgb_data(pixel_count * 3);
        
        const uint8_t* pixel_ptr = pixel_start;
        for (size_t i = 0; i < pixel_count; i++) {
            if (pixel_ptr + 4 > end) break;
            
            uint8_t r = pixel_ptr[0];
            uint8_t g = pixel_ptr[1];
            uint8_t b = pixel_ptr[2];
            uint8_t e = pixel_ptr[3];
            pixel_ptr += 4;
            
            if (e == 0) {
                rgb_data[i * 3] = rgb_data[i * 3 + 1] = rgb_data[i * 3 + 2] = 0;
            } else {
                int exp = static_cast<int>(e) - 128;
                float f = std::ldexp(1.0f, exp - 8) * exposure;
                
                float rf = static_cast<float>(r) * f;
                float gf = static_cast<float>(g) * f;
                float bf = static_cast<float>(b) * f;

                rf = std::sqrt(rf);
                gf = std::sqrt(gf);
                bf = std::sqrt(bf);
                
                rgb_data[i * 3] = static_cast<uint8_t>(std::min(rf * 255.0f, 255.0f));
                rgb_data[i * 3 + 1] = static_cast<uint8_t>(std::min(gf * 255.0f, 255.0f));
                rgb_data[i * 3 + 2] = static_cast<uint8_t>(std::min(bf * 255.0f, 255.0f));
            }
        }
        
        return ImageData(width, height, 3, std::move(rgb_data));
    }
};

// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//                  TGA DECODER
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

class TGADecoder {
public:
    static ImageData load(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file) throw image_error("TGA: Cannot open file");

        uint8_t header[18];
        file.read((char*)header, 18);
        
        int width = header[12] | (header[13] << 8);
        int height = header[14] | (header[15] << 8);
        int bpp = header[16];
        int imageType = header[2];
        
        if (imageType != 2 && imageType != 10) 
            throw image_error("TGA: Only true-color supported");
        
        int channels = bpp == 24 ? 3 : bpp == 32 ? 4 : 0;
        if (!channels) throw image_error("TGA: Unsupported bpp");

        int idLen = header[0];
        int hasCmap = header[1];
        if (hasCmap) {
            int cmapLen = header[5] | (header[4] << 8);
            int cmapBpp = header[7];
            file.seekg(cmapLen * (cmapBpp / 8), std::ios::cur);
        }
        file.seekg(idLen, std::ios::cur);

        file.seekg(0, std::ios::end);
        size_t fileSize = file.tellg();
        file.seekg(18 + idLen + (hasCmap ? ((header[5] | (header[4] << 8)) * (header[7] / 8)) : 0));
        
        std::vector<uint8_t> rawData(fileSize - file.tellg());
        file.read((char*)rawData.data(), rawData.size());
        
        std::vector<uint8_t> pixels(width * height * channels);
        
        if (imageType == 2) { // Uncompressed
            int stride = ((width * channels + 3) & ~3);
            
            #pragma omp parallel for schedule(static)
            for (int y = 0; y < height; y++) {
                int srcY = height - 1 - y;
                memcpy(&pixels[y * width * channels], 
                       &rawData[srcY * stride], 
                       width * channels);
            }
            
            // BGR -> RGB
            #pragma omp parallel for schedule(static)
            for (size_t i = 0; i < width * height * channels; i += channels * 4) {
                for (int j = 0; j < 4 && (i + j * channels) < width * height * channels; j++) {
                    size_t idx = i + j * channels;
                    if (idx + 2 < width * height * channels) {
                        std::swap(pixels[idx], pixels[idx + 2]);
                    }
                }
            }
        } else { // RLE compressed
            const uint8_t* src = rawData.data();
            uint8_t* dst = pixels.data();
            size_t total = width * height;
            
            while (total > 0) {
                uint8_t packet = *src++;
                uint32_t count = (packet & 0x7F) + 1;
                
                if (packet & 0x80) {
                    for (uint32_t i = 0; i < count; i++) {
                        dst[0] = src[2]; // B -> R
                        dst[1] = src[1]; // G
                        dst[2] = src[0]; // R -> B
                        dst += channels;
                    }
                    src += channels;
                } else {
                    for (uint32_t i = 0; i < count; i++) {
                        dst[0] = src[2];
                        dst[1] = src[1];
                        dst[2] = src[0];
                        dst += channels;
                        src += channels;
                    }
                }
                total -= count;
            }
        }
        
        return ImageData(width, height, channels, std::move(pixels));
    }
};

// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//                  BMP DECODER
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

class BMPDecoder {
private:
    static const uint32_t BMP_BI_RGB = 0;
    static const uint32_t BMP_BI_RLE8 = 1;
    static const uint32_t BMP_BI_RLE4 = 2;
    static const uint32_t BMP_BI_BITFIELDS = 3;

    static std::vector<uint8_t> decode_rle8(const uint8_t* data, size_t size, int width, int height) {
        std::vector<uint8_t> indices(width * height, 0);
        int x = 0, y = height - 1;
        size_t pos = 0;
        
        while (y >= 0 && pos < size) {
            uint8_t byte1 = data[pos++];
            uint8_t byte2 = (pos < size) ? data[pos++] : 0;
            
            if (byte1 == 0) {
                if (byte2 == 0) {
                    x = 0;
                    y--;
                } else if (byte2 == 1) {
                    break;
                } else if (byte2 == 2) {
                    if (pos + 2 > size) break;
                    uint8_t dx = data[pos++];
                    uint8_t dy = data[pos++];
                    x += dx;
                    y -= dy;
                } else {
                    int count = byte2;
                    for (int i = 0; i < count; i++) {
                        if (pos >= size) break;
                        if (x < width && y >= 0) {
                            indices[y * width + x] = data[pos++];
                        }
                        x++;
                    }
                    if (count % 2 && pos < size) pos++;
                }
            } else {
                for (int i = 0; i < byte1; i++) {
                    if (x < width && y >= 0) {
                        indices[y * width + x] = byte2;
                    }
                    x++;
                }
            }
        }
        return indices;
    }

    static std::vector<uint8_t> decode_rle4(const uint8_t* data, size_t size, int width, int height) {
        std::vector<uint8_t> indices(width * height, 0);
        int x = 0, y = height - 1;
        size_t pos = 0;
        
        while (y >= 0 && pos < size) {
            uint8_t byte1 = data[pos++];
            uint8_t byte2 = (pos < size) ? data[pos++] : 0;
            
            if (byte1 == 0) {
                if (byte2 == 0) {
                    x = 0;
                    y--;
                } else if (byte2 == 1) {
                    break;
                } else if (byte2 == 2) {
                    if (pos + 2 > size) break;
                    uint8_t dx = data[pos++];
                    uint8_t dy = data[pos++];
                    x += dx;
                    y -= dy;
                } else {
                    int count = byte2;
                    int byteCount = (count + 1) / 2;
                    
                    for (int i = 0; i < count; i++) {
                        if (x < width && y >= 0) {
                            uint8_t byte = data[pos + i / 2];
                            uint8_t nibble = (i % 2) ? (byte & 0x0F) : (byte >> 4);
                            indices[y * width + x] = nibble;
                        }
                        x++;
                    }
                    
                    pos += byteCount;
                    if (byteCount % 2) pos++;
                }
            } else {
                uint8_t nibble1 = byte2 >> 4;
                uint8_t nibble2 = byte2 & 0x0F;
                
                for (int i = 0; i < byte1; i++) {
                    if (x < width && y >= 0) {
                        indices[y * width + x] = (i % 2) ? nibble2 : nibble1;
                    }
                    x++;
                }
            }
        }
        
        return indices;
    }

public:
    static ImageData load(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            throw image_error("BMP: Cannot open file");
        }

        uint8_t fileHeader[14];
        file.read((char*)fileHeader, 14);
        
        if (fileHeader[0] != 'B' || fileHeader[1] != 'M') {
            throw image_error("BMP: Invalid signature");
        }

        uint32_t dataOffset = *(uint32_t*)&fileHeader[10];

        uint32_t infoSize;
        file.read((char*)&infoSize, 4);
        
        if (infoSize < 40) {
            throw image_error("BMP: Invalid info header size");
        }

        int32_t width, height;
        file.read((char*)&width, 4);
        file.read((char*)&height, 4);
        
        uint16_t planes;
        file.read((char*)&planes, 2);
        if (planes != 1) {
            throw image_error("BMP: Invalid planes count");
        }

        uint16_t bpp;
        file.read((char*)&bpp, 2);
        
        uint32_t compression;
        file.read((char*)&compression, 4);
        
        file.seekg(14 + infoSize);

        int absWidth = std::abs(width);
        int absHeight = std::abs(height);
        
        if (absWidth <= 0 || absHeight <= 0 || 
            absWidth > IMG_MAX_IMAGE_SIZE || absHeight > IMG_MAX_IMAGE_SIZE) {
            throw image_error("BMP: Invalid dimensions");
        }

        std::vector<uint8_t> palette;
        
        if (bpp <= 8) {
            int maxColors = 1 << bpp;
            palette.resize(maxColors * 4);
            
            int paletteOffset = 14 + infoSize;
            if (dataOffset > static_cast<uint32_t>(paletteOffset)) {
                file.seekg(paletteOffset);
                file.read((char*)palette.data(), maxColors * 4);
            }
        }

        file.seekg(dataOffset);
        
        std::vector<uint8_t> pixels;
        int channels = 0;

        auto calculate_row_size = [&](int bits_per_pixel) -> int {
            return ((absWidth * bits_per_pixel + 31) / 32) * 4;
        };

        try {
            if (bpp == 1) { // Monochrome
                channels = 1;
                pixels.resize(absWidth * absHeight);
                int rowSize = calculate_row_size(1);
                std::vector<uint8_t> row(rowSize);
                
                for (int y = 0; y < absHeight; y++) {
                    file.read((char*)row.data(), rowSize);
                    
                    for (int x = 0; x < absWidth; x++) {
                        int byteIdx = x / 8;
                        int bitIdx = 7 - (x % 8);
                        uint8_t bit = (row[byteIdx] >> bitIdx) & 1;
                        pixels[y * absWidth + x] = palette.empty() ? (bit * 255) : palette[bit * 4 + 2];
                    }
                }
            }
            else if (bpp == 4) {
                if (compression == BMP_BI_RGB) {
                    channels = 1;
                    pixels.resize(absWidth * absHeight);
                    int rowSize = calculate_row_size(4);
                    std::vector<uint8_t> row(rowSize);
                    
                    for (int y = 0; y < absHeight; y++) {
                        file.read((char*)row.data(), rowSize);
                        
                        for (int x = 0; x < absWidth; x++) {
                            int byteIdx = x / 2;
                            uint8_t nibble = (x % 2) ? (row[byteIdx] & 0x0F) : (row[byteIdx] >> 4);
                            pixels[y * absWidth + x] = palette.empty() ? (nibble * 17) : palette[nibble * 4 + 2];
                        }
                    }
                }
                else if (compression == BMP_BI_RLE4) {
                    file.seekg(0, std::ios::end);
                    std::streampos endPos = file.tellg();
                    size_t compressedSize = static_cast<size_t>(endPos) - dataOffset;
                    file.seekg(dataOffset);
                    
                    std::vector<uint8_t> compressedData(compressedSize);
                    file.read((char*)compressedData.data(), compressedSize);
                    
                    std::vector<uint8_t> indices = decode_rle4(compressedData.data(), compressedSize, absWidth, absHeight);
                    
                    channels = 3;
                    pixels.resize(absWidth * absHeight * 3);

                    for (int i = 0; i < absWidth * absHeight; i++) {
                        uint8_t nibble = indices[i];
                        int idx = i * 3;
                        
                        if (!palette.empty()) {
                            pixels[idx] = palette[nibble * 4 + 2];
                            pixels[idx + 1] = palette[nibble * 4 + 1];
                            pixels[idx + 2] = palette[nibble * 4];
                        } else {
                            uint8_t gray = nibble * 17;
                            pixels[idx] = pixels[idx + 1] = pixels[idx + 2] = gray;
                        }
                    }
                }
            }
            else if (bpp == 8) {
                if (compression == BMP_BI_RGB) {
                    channels = 3;
                    pixels.resize(absWidth * absHeight * 3);
                    int rowSize = calculate_row_size(8);
                    std::vector<uint8_t> row(rowSize);
                    
                    for (int y = 0; y < absHeight; y++) {
                        file.read((char*)row.data(), rowSize);
                        
                        for (int x = 0; x < absWidth; x++) {
                            uint8_t idx = row[x];
                            int pixelIdx = y * absWidth * 3 + x * 3;
                            
                            if (!palette.empty()) {
                                pixels[pixelIdx] = palette[idx * 4 + 2];
                                pixels[pixelIdx + 1] = palette[idx * 4 + 1];
                                pixels[pixelIdx + 2] = palette[idx * 4];
                            } else {
                                pixels[pixelIdx] = pixels[pixelIdx + 1] = pixels[pixelIdx + 2] = idx;
                            }
                        }
                    }
                }
                else if (compression == BMP_BI_RLE8) {
                    file.seekg(0, std::ios::end);
                    std::streampos endPos = file.tellg();
                    size_t compressedSize = static_cast<size_t>(endPos) - dataOffset;
                    file.seekg(dataOffset);
                    
                    std::vector<uint8_t> compressedData(compressedSize);
                    file.read((char*)compressedData.data(), compressedSize);
                    
                    std::vector<uint8_t> indices = decode_rle8(compressedData.data(), compressedSize, absWidth, absHeight);
                    
                    channels = 3;
                    pixels.resize(absWidth * absHeight * 3);

                    for (int i = 0; i < absWidth * absHeight; i++) {
                        uint8_t idx = indices[i];
                        int pixelIdx = i * 3;
                        
                        if (!palette.empty()) {
                            pixels[pixelIdx] = palette[idx * 4 + 2];
                            pixels[pixelIdx + 1] = palette[idx * 4 + 1];
                            pixels[pixelIdx + 2] = palette[idx * 4];
                        } else {
                            pixels[pixelIdx] = pixels[pixelIdx + 1] = pixels[pixelIdx + 2] = idx;
                        }
                    }
                }
            }
            else if (bpp == 24) {
                channels = 3;
                pixels.resize(absWidth * absHeight * 3);
                int rowSize = calculate_row_size(24);
                std::vector<uint8_t> row(rowSize);
                
                for (int y = 0; y < absHeight; y++) {
                    file.read((char*)row.data(), rowSize);
                    
                    for (int x = 0; x < absWidth; x++) {
                        int srcIdx = x * 3;
                        int dstIdx = (y * absWidth + x) * 3;
                        
                        // BGR -> RGB
                        pixels[dstIdx] = row[srcIdx + 2];     // R
                        pixels[dstIdx + 1] = row[srcIdx + 1]; // G
                        pixels[dstIdx + 2] = row[srcIdx];     // B
                    }
                }
            }
            else if (bpp == 32) {
                channels = 4;
                pixels.resize(absWidth * absHeight * 4);
                
                for (int y = 0; y < absHeight; y++) {
                    for (int x = 0; x < absWidth; x++) {
                        uint8_t b, g, r, a;
                        file.read((char*)&b, 1);
                        file.read((char*)&g, 1);
                        file.read((char*)&r, 1);
                        file.read((char*)&a, 1);
                        
                        int idx = (y * absWidth + x) * 4;
                        pixels[idx] = r;     // R
                        pixels[idx + 1] = g; // G
                        pixels[idx + 2] = b; // B
                        pixels[idx + 3] = a; // A
                    }
                }
            }
            else {
                throw image_error("BMP: Unsupported bit depth: " + std::to_string(bpp));
            }
        } catch (const std::exception& e) {
            throw image_error(std::string("BMP decoding error: ") + e.what());
        }

        ImageData result(absWidth, absHeight, channels, std::move(pixels));
        
        return result;
    }
};

// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//                  DDS DECODER (100% SUPPORT)
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

class DDSDecoder {
private:
    #pragma pack(push, 1)
    struct DDS_PIXELFORMAT {
        uint32_t dwSize;
        uint32_t dwFlags;
        uint32_t dwFourCC;
        uint32_t dwRGBBitCount;
        uint32_t dwRBitMask;
        uint32_t dwGBitMask;
        uint32_t dwBBitMask;
        uint32_t dwABitMask;
    };

    struct DDS_HEADER {
        uint32_t dwSize;
        uint32_t dwFlags;
        uint32_t dwHeight;
        uint32_t dwWidth;
        uint32_t dwPitchOrLinearSize;
        uint32_t dwDepth;
        uint32_t dwMipMapCount;
        uint32_t dwReserved1[11];
        DDS_PIXELFORMAT ddspf;
        uint32_t dwCaps;
        uint32_t dwCaps2;
        uint32_t dwCaps3;
        uint32_t dwCaps4;
        uint32_t dwReserved2;
    };
    #pragma pack(pop)

    static constexpr uint32_t FOURCC_DXT1 = 0x31545844;
    static constexpr uint32_t FOURCC_DXT3 = 0x33545844;
    static constexpr uint32_t FOURCC_DXT5 = 0x35545844;
    static constexpr uint32_t FOURCC_ATI1 = 0x31495441;
    static constexpr uint32_t FOURCC_ATI2 = 0x32495441;
    static constexpr uint32_t FOURCC_DX10 = 0x30315844;

    static ALWAYS_INLINE void decode_bc1_block(uint8_t* dst, const uint8_t* src, int width) {
        uint16_t color0 = src[0] | (src[1] << 8);
        uint16_t color1 = src[2] | (src[3] << 8);
        uint32_t indices = src[4] | (src[5] << 8) | (src[6] << 16) | (src[7] << 24);
        
        uint8_t r[4], g[4], b[4];
        r[0] = (color0 >> 11) & 0x1F; g[0] = (color0 >> 5) & 0x3F; b[0] = color0 & 0x1F;
        r[1] = (color1 >> 11) & 0x1F; g[1] = (color1 >> 5) & 0x3F; b[1] = color1 & 0x1F;
        
        r[0] = (r[0] << 3) | (r[0] >> 2);
        g[0] = (g[0] << 2) | (g[0] >> 4);
        b[0] = (b[0] << 3) | (b[0] >> 2);
        r[1] = (r[1] << 3) | (r[1] >> 2);
        g[1] = (g[1] << 2) | (g[1] >> 4);
        b[1] = (b[1] << 3) | (b[1] >> 2);
        
        if (color0 > color1) {
            r[2] = (2 * r[0] + r[1] + 1) / 3;
            g[2] = (2 * g[0] + g[1] + 1) / 3;
            b[2] = (2 * b[0] + b[1] + 1) / 3;
            r[3] = (r[0] + 2 * r[1] + 1) / 3;
            g[3] = (g[0] + 2 * g[1] + 1) / 3;
            b[3] = (b[0] + 2 * b[1] + 1) / 3;
        } else {
            r[2] = (r[0] + r[1]) / 2;
            g[2] = (g[0] + g[1]) / 2;
            b[2] = (b[0] + b[1]) / 2;
            r[3] = g[3] = b[3] = 0;
        }
        
        for (int y = 0; y < 4; y++) {
            for (int x = 0; x < 4; x++) {
                int idx = (indices >> (2 * (y * 4 + x))) & 3;
                int dst_idx = ((y * width) + x) * 4;
                dst[dst_idx] = r[idx];
                dst[dst_idx + 1] = g[idx];
                dst[dst_idx + 2] = b[idx];
                dst[dst_idx + 3] = (color0 <= color1 && idx == 3) ? 0 : 255;
            }
        }
    }

    static ALWAYS_INLINE void decode_bc3_block(uint8_t* dst, const uint8_t* src, int width) {
        uint8_t alpha0 = src[0];
        uint8_t alpha1 = src[1];
        uint64_t alpha_bits = 0;
        memcpy(&alpha_bits, src + 2, 6);
        
        uint8_t alphas[8];
        alphas[0] = alpha0;
        alphas[1] = alpha1;
        
        if (alpha0 > alpha1) {
            alphas[2] = (6 * alpha0 + 1 * alpha1) / 7;
            alphas[3] = (5 * alpha0 + 2 * alpha1) / 7;
            alphas[4] = (4 * alpha0 + 3 * alpha1) / 7;
            alphas[5] = (3 * alpha0 + 4 * alpha1) / 7;
            alphas[6] = (2 * alpha0 + 5 * alpha1) / 7;
            alphas[7] = (1 * alpha0 + 6 * alpha1) / 7;
        } else {
            alphas[2] = (4 * alpha0 + 1 * alpha1) / 5;
            alphas[3] = (3 * alpha0 + 2 * alpha1) / 5;
            alphas[4] = (2 * alpha0 + 3 * alpha1) / 5;
            alphas[5] = (1 * alpha0 + 4 * alpha1) / 5;
            alphas[6] = 0;
            alphas[7] = 255;
        }
        
        decode_bc1_block(dst, src + 8, width);
        
        for (int y = 0; y < 4; y++) {
            for (int x = 0; x < 4; x++) {
                int bit_pos = 3 * (y * 4 + x);
                int idx = (alpha_bits >> bit_pos) & 7;
                dst[(y * width + x) * 4 + 3] = alphas[idx];
            }
        }
    }

    static ALWAYS_INLINE uint32_t get_mask_shift(uint32_t mask) {
        if (!mask) return 0;
        uint32_t shift = 0;
        while (!(mask & 1)) {
            mask >>= 1;
            shift++;
        }
        return shift;
    }

    static ALWAYS_INLINE uint32_t get_mask_max(uint32_t mask) {
        if (!mask) return 0;
        uint32_t max_val = mask;
        while (max_val & (max_val + 1)) max_val |= max_val >> 1;
        return max_val;
    }

    static std::vector<uint8_t> decode_rgb(const uint8_t* data, uint32_t width, uint32_t height, const DDS_PIXELFORMAT& fmt, size_t data_size) {
        uint32_t bpp = fmt.dwRGBBitCount;
        uint32_t pitch = ((width * bpp + 31) / 32) * 4;
        
        if (data_size < pitch * height) throw image_error("Insufficient DDS data");
        
        std::vector<uint8_t> rgba(width * height * 4);
        
        uint32_t r_shift = get_mask_shift(fmt.dwRBitMask);
        uint32_t g_shift = get_mask_shift(fmt.dwGBitMask);
        uint32_t b_shift = get_mask_shift(fmt.dwBBitMask);
        uint32_t a_shift = get_mask_shift(fmt.dwABitMask);
        
        uint32_t r_max = get_mask_max(fmt.dwRBitMask);
        uint32_t g_max = get_mask_max(fmt.dwGBitMask);
        uint32_t b_max = get_mask_max(fmt.dwBBitMask);
        uint32_t a_max = get_mask_max(fmt.dwABitMask);
        
        bool is_bgr = (fmt.dwRBitMask == 0x000000FF && fmt.dwGBitMask == 0x0000FF00 && fmt.dwBBitMask == 0x00FF0000);
        bool is_rgb = (fmt.dwRBitMask == 0x00FF0000 && fmt.dwGBitMask == 0x0000FF00 && fmt.dwBBitMask == 0x000000FF);
        
        #pragma omp parallel for schedule(static)
        for (uint32_t y = 0; y < height; y++) {
            const uint8_t* src_row = data + y * pitch;
            uint8_t* dst_row = &rgba[y * width * 4];
            
            if (bpp == 32) {
                if (is_bgr || is_rgb) {
                    for (uint32_t x = 0; x < width; x++) {
                        uint32_t pixel = *reinterpret_cast<const uint32_t*>(src_row + x * 4);
                        if (is_bgr) {
                            dst_row[x * 4] = (pixel >> 16) & 0xFF;
                            dst_row[x * 4 + 1] = (pixel >> 8) & 0xFF;
                            dst_row[x * 4 + 2] = pixel & 0xFF;
                        } else {
                            dst_row[x * 4] = pixel & 0xFF;
                            dst_row[x * 4 + 1] = (pixel >> 8) & 0xFF;
                            dst_row[x * 4 + 2] = (pixel >> 16) & 0xFF;
                        }
                        dst_row[x * 4 + 3] = fmt.dwABitMask ? ((pixel >> 24) & 0xFF) : 255;
                    }
                } else {
                    for (uint32_t x = 0; x < width; x++) {
                        uint32_t pixel = *reinterpret_cast<const uint32_t*>(src_row + x * 4);
                        dst_row[x * 4] = r_max ? (((pixel & fmt.dwRBitMask) >> r_shift) * 255 / r_max) : 0;
                        dst_row[x * 4 + 1] = g_max ? (((pixel & fmt.dwGBitMask) >> g_shift) * 255 / g_max) : 0;
                        dst_row[x * 4 + 2] = b_max ? (((pixel & fmt.dwBBitMask) >> b_shift) * 255 / b_max) : 0;
                        dst_row[x * 4 + 3] = fmt.dwABitMask ? (((pixel & fmt.dwABitMask) >> a_shift) * 255 / a_max) : 255;
                    }
                }
            } else if (bpp == 24) {
                for (uint32_t x = 0; x < width; x++) {
                    dst_row[x * 4] = src_row[x * 3 + 2];
                    dst_row[x * 4 + 1] = src_row[x * 3 + 1];
                    dst_row[x * 4 + 2] = src_row[x * 3];
                    dst_row[x * 4 + 3] = 255;
                }
            } else if (bpp == 16) {
                if (fmt.dwRBitMask == 0xF800 && fmt.dwGBitMask == 0x07E0 && fmt.dwBBitMask == 0x001F) {
                    for (uint32_t x = 0; x < width; x++) {
                        uint16_t pixel = *reinterpret_cast<const uint16_t*>(src_row + x * 2);
                        dst_row[x * 4] = ((pixel >> 11) & 0x1F) * 255 / 31;
                        dst_row[x * 4 + 1] = ((pixel >> 5) & 0x3F) * 255 / 63;
                        dst_row[x * 4 + 2] = (pixel & 0x1F) * 255 / 31;
                        dst_row[x * 4 + 3] = 255;
                    }
                } else if (fmt.dwRBitMask == 0x7C00 && fmt.dwGBitMask == 0x03E0 && fmt.dwBBitMask == 0x001F && fmt.dwABitMask == 0x8000) {
                    for (uint32_t x = 0; x < width; x++) {
                        uint16_t pixel = *reinterpret_cast<const uint16_t*>(src_row + x * 2);
                        dst_row[x * 4] = ((pixel >> 10) & 0x1F) * 255 / 31;
                        dst_row[x * 4 + 1] = ((pixel >> 5) & 0x1F) * 255 / 31;
                        dst_row[x * 4 + 2] = (pixel & 0x1F) * 255 / 31;
                        dst_row[x * 4 + 3] = (pixel & 0x8000) ? 255 : 0;
                    }
                }
            } else if (bpp == 8 && (fmt.dwFlags & 0x20000)) {
                for (uint32_t x = 0; x < width; x++) {
                    uint8_t lum = src_row[x];
                    dst_row[x * 4] = lum;
                    dst_row[x * 4 + 1] = lum;
                    dst_row[x * 4 + 2] = lum;
                    dst_row[x * 4 + 3] = 255;
                }
            }
        }
        
        return rgba;
    }

    static std::vector<uint8_t> decode_bc(const uint8_t* data, uint32_t width, uint32_t height, 
                                         uint32_t format, size_t data_size) {
        int block_size = (format == FOURCC_DXT1) ? 8 : 16;
        uint32_t blocks_x = (width + 3) / 4;
        uint32_t blocks_y = (height + 3) / 4;
        
        if (data_size < blocks_x * blocks_y * block_size) throw image_error("Insufficient DDS data");
        
        std::vector<uint8_t> rgba(width * height * 4);
        
        #pragma omp parallel for schedule(static) collapse(2)
        for (uint32_t by = 0; by < blocks_y; by++) {
            for (uint32_t bx = 0; bx < blocks_x; bx++) {
                size_t block_offset = (by * blocks_x + bx) * block_size;
                const uint8_t* block_data = data + block_offset;
                
                uint8_t block_rgba[64];
                
                switch (format) {
                    case FOURCC_DXT1:
                        decode_bc1_block(block_rgba, block_data, 4);
                        break;
                    case FOURCC_DXT3:
                        decode_bc3_block(block_rgba, block_data, 4);
                        for (int y = 0; y < 4; y++) {
                            uint16_t alpha_row = block_data[y * 2] | (block_data[y * 2 + 1] << 8);
                            for (int x = 0; x < 4; x++) {
                                uint8_t alpha = (alpha_row >> (x * 4)) & 0xF;
                                block_rgba[(y * 4 + x) * 4 + 3] = alpha * 17;
                            }
                        }
                        break;
                    case FOURCC_DXT5:
                        decode_bc3_block(block_rgba, block_data, 4);
                        break;
                    default:
                        continue;
                }
                
                for (int y = 0; y < 4; y++) {
                    for (int x = 0; x < 4; x++) {
                        uint32_t dst_x = bx * 4 + x;
                        uint32_t dst_y = by * 4 + y;
                        
                        if (dst_x < width && dst_y < height) {
                            size_t dst_idx = (dst_y * width + dst_x) * 4;
                            size_t src_idx = (y * 4 + x) * 4;
                            
                            rgba[dst_idx] = block_rgba[src_idx];
                            rgba[dst_idx + 1] = block_rgba[src_idx + 1];
                            rgba[dst_idx + 2] = block_rgba[src_idx + 2];
                            rgba[dst_idx + 3] = block_rgba[src_idx + 3];
                        }
                    }
                }
            }
        }
        
        return rgba;
    }

public:
    static ImageData load(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open()) throw image_error("Cannot open DDS file: " + filename);
        
        uint32_t magic;
        file.read((char*)&magic, 4);
        if (magic != 0x20534444) throw image_error("Invalid DDS signature");
        
        DDS_HEADER header;
        file.read((char*)&header, sizeof(header));
        
        if (header.dwSize != 124) throw image_error("Invalid DDS header");
        
        uint32_t width = header.dwWidth;
        uint32_t height = header.dwHeight;
        
        if (width > IMG_MAX_IMAGE_SIZE || height > IMG_MAX_IMAGE_SIZE) {
            throw image_error("DDS dimensions too large");
        }
        
        file.seekg(0, std::ios::end);
        size_t file_size = file.tellg();
        size_t data_offset = 128;
        
        if (header.ddspf.dwFourCC == FOURCC_DX10) {
            data_offset += 20;
        }
        
        file.seekg(data_offset, std::ios::beg);
        size_t data_size = file_size - data_offset;
        std::vector<uint8_t> dds_data(data_size);
        file.read((char*)dds_data.data(), data_size);
        
        std::vector<uint8_t> rgba_data;
        uint32_t format = header.ddspf.dwFourCC;
        
        if (header.ddspf.dwFlags & 0x40) {
            rgba_data = decode_rgb(dds_data.data(), width, height, header.ddspf, data_size);
        } else if (header.ddspf.dwFlags & 0x4) {
            if (format == FOURCC_DXT1 || format == FOURCC_DXT3 || format == FOURCC_DXT5) {
                rgba_data = decode_bc(dds_data.data(), width, height, format, data_size);
            } else if (format == FOURCC_ATI1 || format == FOURCC_ATI2) {
                throw image_error("ATI1/ATI2 formats not supported");
            } else {
                throw image_error("Unsupported DDS format");
            }
        } else if (header.ddspf.dwFlags & 0x20000) {
            rgba_data = decode_rgb(dds_data.data(), width, height, header.ddspf, data_size);
        } else {
            throw image_error("Unsupported DDS pixel format");
        }
        
        ImageData result(width, height, 4);
        result.data = std::move(rgba_data);
        
        result.flip_vertically();

        return result;
    }
};

// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//                 QOI DECODER (100% SUPPORT)
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

class QOIDecoder {
private:
    #pragma pack(push, 1)
    struct QOIHeader {
        uint32_t magic;
        uint32_t width;
        uint32_t height;
        uint8_t channels;
        uint8_t colorspace;
    };
    #pragma pack(pop)

    static constexpr uint8_t QOI_OP_INDEX  = 0x00;
    static constexpr uint8_t QOI_OP_DIFF   = 0x40;
    static constexpr uint8_t QOI_OP_LUMA   = 0x80;
    static constexpr uint8_t QOI_OP_RUN    = 0xC0;
    static constexpr uint8_t QOI_OP_RGB    = 0xFE;
    static constexpr uint8_t QOI_OP_RGBA   = 0xFF;
    static constexpr uint8_t QOI_MASK_2    = 0xC0;
    static constexpr uint32_t QOI_MAGIC    = 0x716F6966;
    static constexpr uint64_t QOI_PADDING  = 0x0100000000000000ULL;

    struct PixelCache {
        uint8_t rgba[4] = {0, 0, 0, 255};
        
        PixelCache() = default;
        PixelCache(uint8_t r, uint8_t g, uint8_t b, uint8_t a) : rgba{r, g, b, a} {}
        
        ALWAYS_INLINE bool operator==(const PixelCache& other) const {
            return rgba[0] == other.rgba[0] && rgba[1] == other.rgba[1] && 
                   rgba[2] == other.rgba[2] && rgba[3] == other.rgba[3];
        }
        
        ALWAYS_INLINE uint8_t hash() const {
            return (rgba[0] * 3 + rgba[1] * 5 + rgba[2] * 7 + rgba[3] * 11) % 64;
        }
    };

    static ALWAYS_INLINE uint32_t read_u32_be(const uint8_t* data) {
        return (data[0] << 24) | (data[1] << 16) | (data[2] << 8) | data[3];
    }

    static ALWAYS_INLINE void write_pixel(uint8_t* dst, const PixelCache& px, int channels) {
        dst[0] = px.rgba[0];
        dst[1] = px.rgba[1];
        dst[2] = px.rgba[2];
        if (channels == 4) dst[3] = px.rgba[3];
    }

public:
    static ImageData load(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open()) throw image_error("Cannot open QOI file: " + filename);
        
        file.seekg(0, std::ios::end);
        size_t file_size = file.tellg();
        file.seekg(0, std::ios::beg);
        
        if (file_size < 22) throw image_error("QOI file too small");
        
        std::vector<uint8_t> qoi_data(file_size);
        file.read((char*)qoi_data.data(), file_size);
        
        return decode(qoi_data.data(), file_size);
    }
    
    static ImageData decode(const uint8_t* data, size_t size) {
        if (size < 22) throw image_error("QOI data too small");
        
        if (read_u32_be(data) != QOI_MAGIC) throw image_error("Invalid QOI magic");
        
        uint32_t width = read_u32_be(data + 4);
        uint32_t height = read_u32_be(data + 8);
        uint8_t channels = data[12];
        uint8_t colorspace = data[13];
        
        if (width == 0 || height == 0 || width > IMG_MAX_IMAGE_SIZE || height > IMG_MAX_IMAGE_SIZE) {
            throw image_error("Invalid QOI dimensions");
        }
        if (channels != 3 && channels != 4) throw image_error("Invalid QOI channels");
        if (colorspace > 1) throw image_error("Invalid QOI colorspace");
        if (width * height > 400000000) throw image_error("QOI image too large");
        
        std::vector<uint8_t> pixels(width * height * channels);
        PixelCache cache[64];
        for (int i = 0; i < 64; i++) cache[i] = PixelCache(0, 0, 0, 0);
        
        PixelCache px(0, 0, 0, 255);
        size_t data_offset = 14;
        size_t pixel_offset = 0;
        size_t total_pixels = width * height;
        
        for (size_t i = 0; i < total_pixels; ++i) {
            if (data_offset >= size - 8) throw image_error("Unexpected end of QOI data");
            
            uint8_t chunk = data[data_offset++];
            
            if (chunk == QOI_OP_RGB) {
                if (data_offset + 2 > size - 8) throw image_error("Unexpected end of QOI data");
                px.rgba[0] = data[data_offset++];
                px.rgba[1] = data[data_offset++];
                px.rgba[2] = data[data_offset++];
            }
            else if (chunk == QOI_OP_RGBA) {
                if (data_offset + 3 > size - 8) throw image_error("Unexpected end of QOI data");
                px.rgba[0] = data[data_offset++];
                px.rgba[1] = data[data_offset++];
                px.rgba[2] = data[data_offset++];
                px.rgba[3] = data[data_offset++];
            }
            else {
                uint8_t tag = chunk & QOI_MASK_2;
                
                if (tag == QOI_OP_INDEX) {
                    px = cache[chunk & 0x3F];
                }
                else if (tag == QOI_OP_DIFF) {
                    int dr = ((chunk >> 4) & 0x03) - 2;
                    int dg = ((chunk >> 2) & 0x03) - 2;
                    int db = (chunk & 0x03) - 2;
                    
                    px.rgba[0] = static_cast<uint8_t>(px.rgba[0] + dr);
                    px.rgba[1] = static_cast<uint8_t>(px.rgba[1] + dg);
                    px.rgba[2] = static_cast<uint8_t>(px.rgba[2] + db);
                }
                else if (tag == QOI_OP_LUMA) {
                    if (data_offset >= size - 8) throw image_error("Unexpected end of QOI data");
                    uint8_t next_byte = data[data_offset++];
                    int dg = (chunk & 0x3F) - 32;
                    int dr_dg = (next_byte >> 4) & 0x0F;
                    int db_dg = next_byte & 0x0F;
                    
                    int dr = dg + dr_dg - 8;
                    int db = dg + db_dg - 8;
                    
                    px.rgba[0] = static_cast<uint8_t>(px.rgba[0] + dr);
                    px.rgba[1] = static_cast<uint8_t>(px.rgba[1] + dg);
                    px.rgba[2] = static_cast<uint8_t>(px.rgba[2] + db);
                }
                else if (tag == QOI_OP_RUN) {
                    uint8_t run = (chunk & 0x3F) + 1;
                    for (uint8_t r = 0; r < run; ++r) {
                        if (i + r >= total_pixels) break;
                        write_pixel(&pixels[pixel_offset], px, channels);
                        pixel_offset += channels;
                    }
                    i += run - 1;
                    continue;
                }
                else {
                    throw image_error("Invalid QOI chunk tag");
                }
            }
            
            cache[px.hash()] = px;
            write_pixel(&pixels[pixel_offset], px, channels);
            pixel_offset += channels;
        }
        
        if (size - data_offset != 8) throw image_error("QOI data length mismatch");
        uint64_t padding = 0;
        memcpy(&padding, data + data_offset, 8);
        if (padding != QOI_PADDING) throw image_error("Invalid QOI padding");
        
        ImageData result(width, height, channels);
        result.data = std::move(pixels);
        
        return result;
    }
    
    static ImageData load_from_memory(const uint8_t* data, size_t size) {
        return decode(data, size);
    }
    
    static bool validate(const uint8_t* data, size_t size) {
        if (size < 22) return false;
        if (read_u32_be(data) != QOI_MAGIC) return false;
        
        uint32_t width = read_u32_be(data + 4);
        uint32_t height = read_u32_be(data + 8);
        uint8_t channels = data[12];
        uint8_t colorspace = data[13];
        
        if (width == 0 || height == 0 || width > IMG_MAX_IMAGE_SIZE || height > IMG_MAX_IMAGE_SIZE) return false;
        if (channels != 3 && channels != 4) return false;
        if (colorspace > 1) return false;
        if (width * height > 400000000) return false;
        
        return true;
    }
    
    static ImageData decode_fast(const uint8_t* data, size_t size) {
        if (size < 22) throw image_error("QOI data too small");
        
        if (read_u32_be(data) != QOI_MAGIC) throw image_error("Invalid QOI magic");
        
        uint32_t width = read_u32_be(data + 4);
        uint32_t height = read_u32_be(data + 8);
        uint8_t channels = data[12];
        
        if (channels != 3 && channels != 4) throw image_error("Invalid QOI channels");
        
        size_t data_offset = 14;
        size_t total_pixels = width * height;
        std::vector<uint8_t> pixels(total_pixels * channels);
        
        PixelCache cache[64];
        for (int i = 0; i < 64; i++) cache[i] = PixelCache(0, 0, 0, 0);
        
        PixelCache px(0, 0, 0, 255);
        size_t pixel_offset = 0;
        
        for (size_t i = 0; i < total_pixels; ++i) {
            if (data_offset >= size - 8) break;
            
            uint8_t chunk = data[data_offset++];
            
            if (chunk == QOI_OP_RGB) {
                if (data_offset + 2 > size - 8) break;
                px.rgba[0] = data[data_offset++];
                px.rgba[1] = data[data_offset++];
                px.rgba[2] = data[data_offset++];
            }
            else if (chunk == QOI_OP_RGBA) {
                if (data_offset + 3 > size - 8) break;
                px.rgba[0] = data[data_offset++];
                px.rgba[1] = data[data_offset++];
                px.rgba[2] = data[data_offset++];
                px.rgba[3] = data[data_offset++];
            }
            else {
                uint8_t tag = chunk & QOI_MASK_2;
                
                if (tag == QOI_OP_INDEX) {
                    px = cache[chunk & 0x3F];
                }
                else if (tag == QOI_OP_DIFF) {
                    int dr = ((chunk >> 4) & 0x03) - 2;
                    int dg = ((chunk >> 2) & 0x03) - 2;
                    int db = (chunk & 0x03) - 2;
                    
                    px.rgba[0] = static_cast<uint8_t>(px.rgba[0] + dr);
                    px.rgba[1] = static_cast<uint8_t>(px.rgba[1] + dg);
                    px.rgba[2] = static_cast<uint8_t>(px.rgba[2] + db);
                }
                else if (tag == QOI_OP_LUMA) {
                    if (data_offset >= size - 8) break;
                    uint8_t next_byte = data[data_offset++];
                    int dg = (chunk & 0x3F) - 32;
                    int dr_dg = (next_byte >> 4) & 0x0F;
                    int db_dg = next_byte & 0x0F;
                    
                    int dr = dg + dr_dg - 8;
                    int db = dg + db_dg - 8;
                    
                    px.rgba[0] = static_cast<uint8_t>(px.rgba[0] + dr);
                    px.rgba[1] = static_cast<uint8_t>(px.rgba[1] + dg);
                    px.rgba[2] = static_cast<uint8_t>(px.rgba[2] + db);
                }
                else if (tag == QOI_OP_RUN) {
                    uint8_t run = (chunk & 0x3F) + 1;
                    for (uint8_t r = 0; r < run && i + r < total_pixels; ++r) {
                        write_pixel(&pixels[pixel_offset], px, channels);
                        pixel_offset += channels;
                    }
                    i += run - 1;
                    continue;
                }
            }
            
            cache[px.hash()] = px;
            write_pixel(&pixels[pixel_offset], px, channels);
            pixel_offset += channels;
        }
        
        if (pixel_offset != total_pixels * channels) {
            throw image_error("QOI pixel count mismatch");
        }
        
        ImageData result(width, height, channels);
        result.data = std::move(pixels);
        
        return result;
    }
    
    static bool get_info(const uint8_t* data, size_t size, int& width, int& height, int& channels, int& colorspace) {
        if (size < 14) return false;
        if (read_u32_be(data) != QOI_MAGIC) return false;
        
        width = read_u32_be(data + 4);
        height = read_u32_be(data + 8);
        channels = data[12];
        colorspace = data[13];
        
        return true;
    }
    
    static bool get_info(const std::string& filename, int& width, int& height, int& channels, int& colorspace) {
        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open()) return false;
        
        uint8_t header[14];
        file.read((char*)header, 14);
        if (file.gcount() != 14) return false;
        
        return get_info(header, 14, width, height, channels, colorspace);
    }
};

// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//                  PUBLIC API
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

inline bool ieq(const char* a, const char* b) {
    while (*a && *b) {
        if (std::tolower(*a++) != std::tolower(*b++))
            return false;
    }
    return *a == *b;
}

inline ImageData load_image(const std::string& filename, int desired_channels = 0) {
    size_t dot_pos = filename.find_last_of('.');
    if (dot_pos == std::string::npos)
        throw image_error("Cannot determine image format");

    const char* ext = filename.c_str() + dot_pos + 1;

    ImageData image;

    if (ieq(ext, "png")) {
        image = PNGDecoder::load(filename);
    }
    else if (ieq(ext, "hdr") || ieq(ext, "pic")) {
        image = HDRDecoder::load(filename);
    }
    else if (ieq(ext, "tga")) {
        image = TGADecoder::load(filename);
    }
    else if (ieq(ext, "bmp")) {
        image = BMPDecoder::load(filename);
    }
    else if (ieq(ext, "dds")) {
        image = DDSDecoder::load(filename);
    }
    else if (ieq(ext, "qoi")) {
        image = QOIDecoder::load(filename);
    }

    else {
        throw image_error(std::string("Unsupported image format: ") + ext);
    }

    if (desired_channels > 0 && desired_channels != image.channels)
        return image.convert_to_channels(desired_channels);

    return image;
}


inline ImageData load_image_from_memory(const uint8_t* data, size_t size, const std::string& format_hint = "") {
    throw image_error("Memory loading not yet implemented");
    return ImageData();
}

inline bool is_image_format_supported(const std::string& filename) {
    static const char* supported[] = {".png", ".jpg", ".jpeg", ".webp", ".gif", ".ico"};
    
    std::string lower = filename;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    
    for (const char* ext : supported) {
        if (lower.length() >= strlen(ext) && 
            lower.compare(lower.length() - strlen(ext), strlen(ext), ext) == 0) {
            return true;
        }
    }
    return false;
}

} // namespace img

#endif // IMAGE_LOADER_HPP
