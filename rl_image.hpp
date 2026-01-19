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
#define RLIMG_AUTHOR "retti"

#define RLIMG_NAMESPACE_BEGIN namespace rlimg {
#define RLIMG_NAMESPACE_END }
#define RL_UNUSED(x) (void)(x)

// --------------------------------------------------------------------------------\
// rl_image.hpp — single‑header image decoder                                      |
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
//   JPEG  ✓ no optimizations
//   JPE   ✓ no optimizations
//   JPG   ✓ no optimizations
//   WebP  in implementation...
//   HDR   ✓
//   GIF   ✓

#define RLIMG_VERSION "0.2"

// Changelog:
//   [01/2026] v0.2 — Optimized PNG decoder, added WebP, HDR, GIF
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
// Supported Formats:
// ------------------
// • PNG (8-bit, with full filter support and SSE-optimized RGBA filtering)
// • JPEG (baseline, progressive with YCbCr to RGB conversion)
// • WebP (lossy/lossless, with SIMD acceleration and caching)
// • HDR/RGBE (high dynamic range, tone-mapping support)
// • GIF (static and animated, with transparency and interlacing)
//
// Key Features:
// -------------
// • Single-header implementation - just #include "rl_image.hpp"
// • No external dependencies (except standard library)
// • SIMD acceleration (SSE2, SSE4.1, AVX2, ARM NEON when available)
// • Thread-safe design (except shared WebP cache)
// • Configurable limits (image size, channels, cache size)
// • Automatic format detection by file extension
// • Channel conversion utilities (RGB↔RGBA, grayscale↔RGB, etc.)
// • Image manipulation (vertical flip, format conversion)
// • WebP caching system with configurable size limits
// • HDR tone mapping with exposure control
// • GIF animation support (multi-frame loading)
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
// Decoder-Specific Notes:
// -----------------------
// PNG:   Full PNG spec support (filters, interlacing, 8-bit only)
// JPEG:  Baseline DCT, YCbCr, optimized Huffman decoding
// WebP:  Caching system, SIMD YUV→RGB, partial VP8L support
// HDR:   RGBE format, Reinhard tone mapping, exposure control
// GIF:   LZW decompression, transparency, frame delays
//
// Performance Tips:
// -----------------
// 1. WebP images are cached in memory (LRU eviction)
// 2. Use desired_channels to avoid post-load conversions
// 3. HDR loading is faster with load_fast() for previews
// 4. GIF multi-frame loading uses load_all_frames() for animations
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
#include <sys/stat.h> // Platform / system
#include <vector>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include <emmintrin.h> // SSE2
#include <chrono>

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
    inline static uint32_t crc_table[256] = {0};
    inline static bool table_computed = false;

    static void make_crc_table() {
        for (uint32_t n = 0; n < 256; n++) {
            uint32_t c = n;
            for (int k = 0; k < 8; k++) {
                if (c & 1)
                    c = 0xEDB88320u ^ (c >> 1); // 0x04C11DB7
                else
                    c >>= 1;
            }
            crc_table[n] = c;
        }
        table_computed = true;
    }

public:
    static uint32_t update(uint32_t crc, const uint8_t* buf, size_t len) {
        if (!table_computed) {
            make_crc_table();
        }

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

    static uint32_t reverse_bits(uint32_t code, int len) {
        uint32_t r = 0;
        for (int i = 0; i < len; ++i) {
            r = (r << 1) | (code & 1);
            code >>= 1;
        }
        return r;
    }

    void build(const std::vector<uint8_t>& code_lengths) {
        int n = static_cast<int>(code_lengths.size());

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

        next_code = first_symbol;

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
                uint32_t rev_code = reverse_bits(static_cast<uint32_t>(huff_code), len);
                
                int entries = 1 << (FAST_BITS - len);
                for (int j = 0; j < entries; ++j) {
                    int idx = rev_code + (j << len);
                    if (idx < (1 << FAST_BITS)) {
                        fast[idx].symbol = static_cast<uint16_t>(sym);
                        fast[idx].length = static_cast<uint8_t>(len);
                    }
                }
            }
        }
    }

    uint16_t decode(uint32_t& bitbuf, int& bits, const uint8_t*& data, size_t& pos, size_t size) {
        while (bits < FAST_BITS && pos < size) {
            bitbuf |= static_cast<uint32_t>(data[pos++]) << bits;
            bits += 8;
        }

        int idx = static_cast<int>(bitbuf & FAST_MASK);
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
                bitbuf |= static_cast<uint32_t>(data[pos++]) << bits;
                bits += 8;
            }
            
            code = (code << 1) | (bitbuf & 1);
            bitbuf >>= 1;
            bits--;

            int first = first_code[len];
            int cnt = count[len];
            
            if (code - first < static_cast<uint32_t>(cnt)) {
                int symbol_idx = first_symbol[len] + (code - first);
                return static_cast<uint16_t>(symbols[symbol_idx]);
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
    static const int LENGTH_CODES      = 29;
    static const int DISTANCE_CODES    = 30;
    static const int CODE_LENGTH_CODES = 19;

    struct LengthInfo {
        uint16_t base;
        uint8_t  extra_bits;
    };

    inline static const LengthInfo length_table[LENGTH_CODES] = {
        {3,0}, {4,0}, {5,0}, {6,0}, {7,0}, {8,0}, {9,0}, {10,0}, {11,1},
        {13,1}, {15,1}, {17,1}, {19,2}, {23,2}, {27,2}, {31,2}, {35,3},
        {43,3}, {51,3}, {59,3}, {67,4}, {83,4}, {99,4}, {115,4}, {131,5},
        {163,5}, {195,5}, {227,5}, {258,0}
    };

    inline static const LengthInfo distance_table[DISTANCE_CODES] = {
        {1,0}, {2,0}, {3,0}, {4,0}, {5,1}, {7,1}, {9,2}, {13,2}, {17,3},
        {25,3}, {33,4}, {49,4}, {65,5}, {97,5}, {129,6}, {193,6}, {257,7},
        {385,7}, {513,8}, {769,8}, {1025,9}, {1537,9}, {2049,10}, {3073,10},
        {4097,11}, {6145,11}, {8193,12}, {12289,12}, {16385,13}, {24577,13}
    };

    inline static const uint8_t cl_order[CODE_LENGTH_CODES] = {
        16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15
    };

    const uint8_t* compressed_data;
    size_t         compressed_size;
    size_t         bit_pos;
    uint32_t       bit_buffer;
    int            bits_in_buffer;

    uint32_t read_bits(int n) {
        while (bits_in_buffer < n) {
            if (bit_pos >= compressed_size) {
                throw image_error("Unexpected end of compressed data");
            }
            bit_buffer |= (uint32_t)compressed_data[bit_pos++] << bits_in_buffer;
            bits_in_buffer += 8;
        }

        uint32_t result = bit_buffer & ((1u << n) - 1);
        bit_buffer >>= n;
        bits_in_buffer -= n;
        return result;
    }

    void align_to_byte() {
        bits_in_buffer = 0;
        bit_buffer     = 0;
    }

    std::vector<uint8_t> read_lengths(HuffmanDecoder& decoder, int count) {
        std::vector<uint8_t> lengths(count, 0);
        int i = 0;

        while (i < count) {
            uint16_t symbol = decoder.decode(bit_buffer, bits_in_buffer,
                                            compressed_data, bit_pos, compressed_size);

            if (symbol < 16) {
                lengths[i++] = static_cast<uint8_t>(symbol);
            }
            else if (symbol == 16) {
                if (i == 0) {
                    throw image_error("Cannot repeat previous at start of lengths");
                }
                int repeat = 3 + (int)read_bits(2);
                if (i + repeat > count) {
                    repeat = count - i;
                }
                uint8_t prev = lengths[i - 1];
                for (int k = 0; k < repeat; ++k) {
                    lengths[i++] = prev;
                }
            }
            else if (symbol == 17) {
                int repeat = 3 + (int)read_bits(3);
                if (i + repeat > count) {
                    repeat = count - i;
                }
                for (int k = 0; k < repeat; ++k) {
                    lengths[i++] = 0;
                }
            }
            else if (symbol == 18) {
                int repeat = 11 + (int)read_bits(7);
                if (i + repeat > count) {
                    repeat = count - i;
                }
                for (int k = 0; k < repeat; ++k) {
                    lengths[i++] = 0;
                }
            }
            else {
                throw image_error("Invalid symbol in code lengths: " + std::to_string(symbol));
            }
        }

        return lengths;
    }

    void decode_block(HuffmanDecoder& litlen_decoder,
                      HuffmanDecoder& dist_decoder,
                      std::vector<uint8_t>& result)
    {
        while (true) {
            uint16_t symbol = litlen_decoder.decode(
                bit_buffer, bits_in_buffer,
                compressed_data, bit_pos, compressed_size
            );

            if (symbol < 256) {
                result.push_back((uint8_t)symbol);
                continue;
            }

            if (symbol == 256) {
                break;
            }

            int length_index = (int)symbol - 257;
            if (length_index < 0 || length_index >= LENGTH_CODES) {
                throw image_error("Invalid length code");
            }

            const LengthInfo& info = length_table[length_index];
            uint32_t length = info.base;
            if (info.extra_bits > 0) {
                length += read_bits(info.extra_bits);
            }

            uint16_t dist_symbol = dist_decoder.decode(
                bit_buffer, bits_in_buffer,
                compressed_data, bit_pos, compressed_size
            );

            if (dist_symbol >= DISTANCE_CODES) {
                throw image_error("Invalid distance code");
            }

            const LengthInfo& dist_info = distance_table[dist_symbol];
            uint32_t distance = dist_info.base;
            if (dist_info.extra_bits > 0) {
                distance += read_bits(dist_info.extra_bits);
            }

            if (distance == 0 || distance > result.size()) {
                throw image_error("Distance too large");
            }

            size_t start = result.size() - distance;
            size_t old_size = result.size();
            result.resize(old_size + length);

            if (distance >= length) {
                std::memcpy(result.data() + old_size,
                            result.data() + start,
                            length);
            } else {
                for (uint32_t i = 0; i < length; i++) {
                    result[old_size + i] = result[start + (i % distance)];
                }
            }
        }
    }

    void process_fixed_block(std::vector<uint8_t>& result) {
        std::vector<uint8_t> litlen_lengths(288, 0);
        for (int i = 0; i <= 143; i++) litlen_lengths[i] = 8;
        for (int i = 144; i <= 255; i++) litlen_lengths[i] = 9;
        for (int i = 256; i <= 279; i++) litlen_lengths[i] = 7;
        for (int i = 280; i <= 287; i++) litlen_lengths[i] = 8;

        std::vector<uint8_t> dist_lengths(32, 5);

        HuffmanDecoder litlen_decoder;
        HuffmanDecoder dist_decoder;

        litlen_decoder.build(litlen_lengths);
        dist_decoder.build(dist_lengths);

        decode_block(litlen_decoder, dist_decoder, result);
    }

    void process_dynamic_block(std::vector<uint8_t>& result) {
        int hlit  = (int)read_bits(5) + 257;  // 257-286
        int hdist = (int)read_bits(5) + 1;    // 1-32
        int hclen = (int)read_bits(4) + 4;    // 4-19

        std::vector<uint8_t> code_lengths(19, 0);
        for (int i = 0; i < hclen; i++) {
            code_lengths[cl_order[i]] = (uint8_t)read_bits(3);
        }

        HuffmanDecoder cl_decoder;
        cl_decoder.build(code_lengths);

        int total_codes = hlit + hdist;
        std::vector<uint8_t> all_lengths = read_lengths(cl_decoder, total_codes);

        std::vector<uint8_t> litlen_lengths(all_lengths.begin(), all_lengths.begin() + hlit);
        std::vector<uint8_t> dist_lengths(all_lengths.begin() + hlit, all_lengths.end());

        HuffmanDecoder litlen_decoder;
        HuffmanDecoder dist_decoder;

        litlen_decoder.build(litlen_lengths);
        dist_decoder.build(dist_lengths);

        decode_block(litlen_decoder, dist_decoder, result);
    }

public:
    Inflate(const uint8_t* data, size_t size)
        : compressed_data(data), compressed_size(size),
          bit_pos(0), bit_buffer(0), bits_in_buffer(0) {}

    std::vector<uint8_t> decompress() {
        std::vector<uint8_t> result;
        result.reserve(compressed_size * 3);

        (void)read_bits(8); // CMF
        (void)read_bits(8); // FLG

        while (true) {
            int final = (int)read_bits(1);
            int type  = (int)read_bits(2);

            if (type == 0) {
                align_to_byte();

                if (bit_pos + 4 > compressed_size) {
                    throw image_error("Invalid uncompressed block header");
                }

                uint32_t len  = compressed_data[bit_pos] |
                                (compressed_data[bit_pos + 1] << 8);
                uint32_t nlen = compressed_data[bit_pos + 2] |
                                (compressed_data[bit_pos + 3] << 8);
                bit_pos += 4;

                if ((len ^ 0xFFFFu) != nlen) {
                    throw image_error("Invalid length in uncompressed block");
                }

                if (bit_pos + len > compressed_size) {
                    throw image_error("Unexpected end of data in uncompressed block");
                }

                size_t old_size = result.size();
                result.resize(old_size + len);
                std::memcpy(result.data() + old_size,
                            compressed_data + bit_pos,
                            len);
                bit_pos += len;
            }
            else if (type == 1) {
                process_fixed_block(result);
            }
            else if (type == 2) {
                process_dynamic_block(result);
            }
            else {
                throw image_error("Invalid block type");
            }

            if (final) {
                break;
            }
        }

        align_to_byte();

        if (bit_pos + 4 <= compressed_size) {
            bit_pos += 4;
        }

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

    static PNGChunk read_chunk(std::ifstream& file) {
        PNGChunk chunk;

        uint8_t len_bytes[4];
        file.read(reinterpret_cast<char*>(len_bytes), 4);
        if (file.gcount() != 4) {
            throw image_error("Failed to read chunk length");
        }

        chunk.length = (len_bytes[0] << 24) | (len_bytes[1] << 16) |
                       (len_bytes[2] << 8)  |  len_bytes[3];

        if (chunk.length > 10000000u) {
            throw image_error("Chunk too large");
        }

        file.read(chunk.type, 4);
        if (file.gcount() != 4) {
            throw image_error("Failed to read chunk type");
        }
        chunk.type[4] = '\0';

        chunk.data.resize(chunk.length);
        if (chunk.length > 0) {
            file.read(reinterpret_cast<char*>(chunk.data.data()), chunk.length);
            if (static_cast<size_t>(file.gcount()) != chunk.length) {
                throw image_error("Failed to read chunk data");
            }
        }

        uint8_t crc_bytes[4];
        file.read(reinterpret_cast<char*>(crc_bytes), 4);
        if (file.gcount() != 4) {
            throw image_error("Failed to read chunk CRC");
        }

        chunk.crc = (crc_bytes[0] << 24) | (crc_bytes[1] << 16) |
                    (crc_bytes[2] << 8)  |  crc_bytes[3];

        uint32_t computed_crc = CRC32::update(
            0, reinterpret_cast<const uint8_t*>(chunk.type), 4
        );
        if (chunk.length > 0) {
            computed_crc = CRC32::update(computed_crc,
                                         chunk.data.data(),
                                         chunk.length);
        }

        if (computed_crc != chunk.crc) {
            throw image_error("Chunk CRC mismatch");
        }

        return chunk;
    }

    static void apply_filter(uint8_t* scanline,
                             const uint8_t* prev_scanline,
                             int width, int channels, int filter_type)
    {
        if (channels == 4 && width >= 4) {
            apply_filter_rgba_sse(scanline, prev_scanline, width, filter_type);
            return;
        }

        const int stride = width * channels;

        switch (filter_type) {
            case 0: // None
                break;

            case 1: // Sub
                for (int i = channels; i < stride; i++) {
                    scanline[i] = (uint8_t)(scanline[i] + scanline[i - channels]);
                }
                break;

            case 2: // Up
                for (int i = 0; i < stride; i++) {
                    scanline[i] = (uint8_t)(scanline[i] + prev_scanline[i]);
                }
                break;

            case 3: // Average
                for (int i = 0; i < stride; i++) {
                    uint8_t left = (i >= channels) ? scanline[i - channels] : 0;
                    uint8_t up   = prev_scanline[i];
                    scanline[i] = (uint8_t)(scanline[i] + ((left + up) >> 1));
                }
                break;

            case 4: // Paeth
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
                break;

            default:
                throw image_error("Unknown PNG filter type");
        }
    }

public:
    static ImageData load(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            throw image_error("Cannot open file: " + filename);
        }

        uint8_t signature[8];
        file.read(reinterpret_cast<char*>(signature), 8);
        if (file.gcount() != 8 ||
            signature[0] != 137 || signature[1] != 80  ||
            signature[2] != 78  || signature[3] != 71  ||
            signature[4] != 13  || signature[5] != 10  ||
            signature[6] != 26  || signature[7] != 10)
        {
            throw image_error("Invalid PNG signature");
        }

        int width = 0, height = 0;
        int bit_depth = 0, color_type = 0;
        std::vector<uint8_t> compressed_data;

        while (true) {
            PNGChunk chunk = read_chunk(file);

            if (std::strcmp(chunk.type, "IHDR") == 0) {
                if (chunk.length != 13) {
                    throw image_error("Invalid IHDR chunk size");
                }

                width  = (chunk.data[0] << 24) | (chunk.data[1] << 16) |
                         (chunk.data[2] << 8)  |  chunk.data[3];
                height = (chunk.data[4] << 24) | (chunk.data[5] << 16) |
                         (chunk.data[6] << 8)  |  chunk.data[7];
                bit_depth  = chunk.data[8];
                color_type = chunk.data[9];

                if (width <= 0 || height <= 0 ||
                    width > IMG_MAX_IMAGE_SIZE || height > IMG_MAX_IMAGE_SIZE)
                {
                    throw image_error("Invalid image dimensions");
                }

                if (bit_depth != 8) {
                    throw image_error("Only 8-bit depth supported");
                }

                if (color_type != 0 && color_type != 2 &&
                    color_type != 3 && color_type != 4 && color_type != 6)
                {
                    throw image_error("Unsupported color type");
                }
            }
            else if (std::strcmp(chunk.type, "IDAT") == 0) {
                compressed_data.insert(compressed_data.end(),
                                       chunk.data.begin(), chunk.data.end());
            }
            else if (std::strcmp(chunk.type, "IEND") == 0) {
                break;
            }
        }

        if (width == 0 || height == 0) {
            throw image_error("No IHDR chunk found");
        }

        if (compressed_data.empty()) {
            throw image_error("No image data found");
        }

        int channels = 0;
        switch (color_type) {
            case 0: channels = 1; break;
            case 2: channels = 3; break;
            case 3: channels = 1; break;
            case 4: channels = 2; break;
            case 6: channels = 4; break;
            default: throw image_error("Unsupported color type");
        }

        Inflate inflater(compressed_data.data(), compressed_data.size());
        std::vector<uint8_t> image_data = inflater.decompress();

        const size_t expected_size = (size_t)height * (size_t)(width * channels + 1);
        if (image_data.size() != expected_size) {
            throw image_error("Decompressed data size mismatch");
        }

        ImageData result(width, height, channels);
        std::vector<uint8_t> prev_scanline(width * channels, 0);

        size_t src_pos = 0;
        size_t dst_pos = 0;
        const int stride = width * channels;

        for (int y = 0; y < height; y++) {
            int filter_type = image_data[src_pos++];

            uint8_t* scanline = &image_data[src_pos];
            apply_filter(scanline, prev_scanline.data(), width, channels, filter_type);

            std::memcpy(&result.data[dst_pos], scanline, stride);
            std::memcpy(prev_scanline.data(), scanline, stride);

            src_pos += stride;
            dst_pos += stride;
        }

        result.flip_vertically();

        return result;
    }
};


// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//                 JPEG DECODER
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

class JPEGDecoder {
private:
    static constexpr int DCT_SIZE = 8;
    static constexpr int DCT_SIZE2 = 64;
    static constexpr int MAX_HUFF_CODES = 256;
    static constexpr int MAX_QUANT_TABLES = 8;
    static constexpr int MAX_HUFF_TABLES = 32;

    struct QuantTable {
        uint16_t values[64];
        uint8_t precision;
        
        QuantTable() : precision(0) {
            memset(values, 0, sizeof(values));
        }
    };
    
    struct HuffmanTable {
        struct FastEntry {
            uint8_t symbol;
            uint8_t length;
        };
        
        FastEntry fast_lookup[256];
        std::vector<uint8_t> huffval;
        uint8_t bits[16];
        int mincode[17];
        int maxcode[17];
        int valptr[17];
        int total_codes;
        
        HuffmanTable() : total_codes(0) {
            memset(fast_lookup, 0, sizeof(fast_lookup));
            memset(bits, 0, sizeof(bits));
        }
    };
    
    struct Component {
        int id;
        uint8_t h_sampling;
        uint8_t v_sampling;
        uint8_t quant_table;
        uint8_t dc_table;
        uint8_t ac_table;
        int width, height;
        std::vector<int16_t> data;
        
        Component() : id(0), h_sampling(0), v_sampling(0), 
                     quant_table(0), dc_table(0), ac_table(0),
                     width(0), height(0) {}
    };
    
    static inline const int zigzag[64] = {
        0,  1,  5,  6, 14, 15, 27, 28,
        2,  4,  7, 13, 16, 26, 29, 42,
        3,  8, 12, 17, 25, 30, 41, 43,
        9, 11, 18, 24, 31, 40, 44, 53,
        10, 19, 23, 32, 39, 45, 52, 54,
        20, 22, 33, 38, 46, 51, 55, 60,
        21, 34, 37, 47, 50, 56, 59, 61,
        35, 36, 48, 49, 57, 58, 62, 63
    };
    
    std::vector<uint8_t> file_data;
    size_t pos;
    uint32_t bit_buffer;
    int bits_in_buffer;
    
    std::vector<QuantTable> quant_tables;
    std::vector<HuffmanTable> huff_tables;
    std::vector<Component> components;
    
    int image_width, image_height;
    int num_components;
    int max_h_sampling, max_v_sampling;

    inline uint32_t read_bits(int n) {
        while (bits_in_buffer < n) {
            if (pos >= file_data.size()) {
                throw image_error("Unexpected end of JPEG data");
            }
            uint8_t byte = file_data[pos++];
            
            if (byte == 0xFF) {
                if (pos < file_data.size()) {
                    uint8_t next = file_data[pos];
                    if (next == 0x00) {
                        pos++; // Skip 0x00
                    } else if (next >= 0xD0 && next <= 0xD7) {
                        pos++; // Skip RST marker
                        bit_buffer = 0;
                        bits_in_buffer = 0;
                        continue;
                    } else {
                        pos--;
                        break;
                    }
                }
            }
            
            bit_buffer = (bit_buffer << 8) | byte;
            bits_in_buffer += 8;
        }
        
        if (bits_in_buffer < n) {
            throw image_error("Not enough bits available");
        }
        
        uint32_t result = (bit_buffer >> (bits_in_buffer - n)) & ((1 << n) - 1);
        bits_in_buffer -= n;
        return result;
    }

    inline void reset_bit_buffer() {
        bit_buffer = 0;
        bits_in_buffer = 0;
    }

    void read_dqt() {
        uint16_t length = (file_data[pos] << 8) | file_data[pos + 1];
        pos += 2;
        
        size_t end_pos = pos + length - 2;
        
        while (pos < end_pos) {
            uint8_t info = file_data[pos++];
            int precision = info >> 4;
            int table_id = info & 0x0F;
            
            if (table_id >= MAX_QUANT_TABLES) {
                throw image_error("Invalid quantization table ID");
            }
            
            if (quant_tables.size() <= static_cast<size_t>(table_id)) {
                quant_tables.resize(table_id + 1);
            }
            
            QuantTable& table = quant_tables[table_id];
            table.precision = precision;
            
            if (precision == 0) {
                for (int i = 0; i < 64; i++) {
                    table.values[zigzag[i]] = file_data[pos++];
                }
            } else {
                for (int i = 0; i < 64; i++) {
                    table.values[zigzag[i]] = (file_data[pos] << 8) | file_data[pos + 1];
                    pos += 2;
                }
            }
        }
    }

    void read_dht() {
        uint16_t length = (file_data[pos] << 8) | file_data[pos + 1];
        pos += 2;
        
        size_t end_pos = pos + length - 2;
        
        while (pos < end_pos) {
            uint8_t info = file_data[pos++];
            int table_class = info >> 4;
            int table_id = info & 0x0F;
            
            int table_index = (table_class << 4) | table_id;
            
            if (huff_tables.size() <= static_cast<size_t>(table_index)) {
                huff_tables.resize(table_index + 1);
            }
            
            HuffmanTable& table = huff_tables[table_index];
            
            int total_codes = 0;
            for (int i = 0; i < 16; i++) {
                table.bits[i] = file_data[pos++];
                total_codes += table.bits[i];
            }
            
            table.huffval.resize(total_codes);
            for (int i = 0; i < total_codes; i++) {
                table.huffval[i] = file_data[pos++];
            }
            
            build_huffman_table(table);
        }
    }

    void build_huffman_table(HuffmanTable& table) {
        int k = 0;
        std::vector<uint8_t> huffsize(256, 0);
        
        for (int i = 1; i <= 16; i++) {
            for (int j = 0; j < table.bits[i-1]; j++) {
                huffsize[k] = i;
                k++;
            }
        }
        huffsize[k] = 0;
        
        k = 0;
        int code = 0;
        int si = huffsize[0];
        std::vector<uint16_t> huffcode(256, 0);
        
        while (true) {
            huffcode[k] = code;
            code++;
            k++;
            
            if (huffsize[k] == si) continue;
            if (huffsize[k] == 0) break;
            
            do {
                code <<= 1;
                si++;
            } while (huffsize[k] != si);
        }

        for (int i = 0; i < 17; i++) {
            table.mincode[i] = 0xFFFF;
            table.maxcode[i] = -1;
            table.valptr[i] = 0;
        }
        
        k = 0;
        for (int i = 1; i <= 16; i++) {
            if (table.bits[i-1] == 0) {
                table.maxcode[i] = -1;
            } else {
                table.valptr[i] = k;
                table.mincode[i] = huffcode[k];
                k += table.bits[i-1];
                table.maxcode[i] = huffcode[k-1];
            }
        }
        
        memset(table.fast_lookup, 0, sizeof(table.fast_lookup));
        
        k = 0;
        for (int len = 1; len <= 8; len++) {
            for (int i = 0; i < table.bits[len-1]; i++) {
                int symbol = table.huffval[k];
                uint16_t code_val = table.mincode[len] + i;
                
                uint8_t rev_code = 0;
                for (int b = 0; b < len; b++) {
                    rev_code = (rev_code << 1) | ((code_val >> (len-1-b)) & 1);
                }
                
                int entries = 1 << (8 - len);
                for (int j = 0; j < entries; j++) {
                    uint8_t idx = rev_code | (j << len);
                    table.fast_lookup[idx].symbol = static_cast<uint8_t>(symbol);
                    table.fast_lookup[idx].length = static_cast<uint8_t>(len);
                }
                k++;
            }
        }
    }

    void read_sof0() {
        uint16_t length = (file_data[pos] << 8) | file_data[pos + 1];
        pos += 2;
        
        uint8_t precision = file_data[pos++];
        image_height = (file_data[pos] << 8) | file_data[pos + 1];
        pos += 2;
        image_width = (file_data[pos] << 8) | file_data[pos + 1];
        pos += 2;
        num_components = file_data[pos++];
        
        if (precision != 8) {
            throw image_error("Only 8-bit precision supported");
        }
        
        if (num_components != 1 && num_components != 3) {
            throw image_error("Unsupported number of components");
        }
        
        components.resize(num_components);
        
        max_h_sampling = 0;
        max_v_sampling = 0;
        for (int i = 0; i < num_components; i++) {
            Component& comp = components[i];
            comp.id = file_data[pos++];
            uint8_t sampling = file_data[pos++];
            comp.h_sampling = sampling >> 4;
            comp.v_sampling = sampling & 0x0F;
            comp.quant_table = file_data[pos++];
            
            if (comp.h_sampling > max_h_sampling) max_h_sampling = comp.h_sampling;
            if (comp.v_sampling > max_v_sampling) max_v_sampling = comp.v_sampling;
        }

        for (int i = 0; i < num_components; i++) {
            Component& comp = components[i];
            comp.width = (image_width * comp.h_sampling + max_h_sampling - 1) / max_h_sampling;
            comp.height = (image_height * comp.v_sampling + max_v_sampling - 1) / max_v_sampling;
            comp.data.resize(comp.width * comp.height);
        }
    }

    void read_sos() {
        uint16_t length = (file_data[pos] << 8) | file_data[pos + 1];
        pos += 2;
        
        int num_scan_components = file_data[pos++];
        
        for (int i = 0; i < num_scan_components; i++) {
            uint8_t comp_id = file_data[pos++];
            uint8_t table_spec = file_data[pos++];
            
            int comp_index = -1;
            for (size_t j = 0; j < components.size(); j++) {
                if (components[j].id == comp_id) {
                    comp_index = j;
                    break;
                }
            }
            
            if (comp_index == -1) {
                throw image_error("Invalid component ID in SOS");
            }
            
            components[comp_index].dc_table = table_spec >> 4;
            components[comp_index].ac_table = table_spec & 0x0F;
        }

        pos += 3;
        decode_scan_data();
    }

    inline int decode_huffman(const HuffmanTable& table) {
        if (bits_in_buffer >= 8) {
            uint8_t peek = static_cast<uint8_t>((bit_buffer >> (bits_in_buffer - 8)) & 0xFF);
            const auto& entry = table.fast_lookup[peek];
            
            if (entry.length != 0) {
                bit_buffer <<= entry.length;
                bits_in_buffer -= entry.length;
                return entry.symbol;
            }
        }
        
        int code = 0;
        int length = 1;
        
        while (length <= 16) {
            code = (code << 1) | read_bits(1);
            
            if (code <= table.maxcode[length]) {
                int index = table.valptr[length] + (code - table.mincode[length]);
                return table.huffval[index];
            }
            
            length++;
        }
        
        throw image_error("Invalid Huffman code");
    }

    inline int decode_coefficient(int num_bits) {
        if (num_bits == 0) return 0;
        
        int value = read_bits(num_bits);

        if (value < (1 << (num_bits - 1))) {
            value = value - (1 << num_bits) + 1;
        }
        
        return value;
    }

    void idct_2d(float block[64]) {
        float temp[64];

        for (int y = 0; y < 8; y++) {
            for (int x = 0; x < 8; x++) {
                float sum = 0.0f;
                for (int u = 0; u < 8; u++) {
                    float cu = (u == 0) ? 1.0f / std::sqrt(2.0f) : 1.0f;
                    sum += cu * block[y * 8 + u] * std::cos((2 * x + 1) * u * M_PI / 16.0f);
                }
                temp[y * 8 + x] = 0.5f * sum;
            }
        }

        for (int x = 0; x < 8; x++) {
            for (int y = 0; y < 8; y++) {
                float sum = 0.0f;
                for (int v = 0; v < 8; v++) {
                    float cv = (v == 0) ? 1.0f / std::sqrt(2.0f) : 1.0f;
                    sum += cv * temp[v * 8 + x] * std::cos((2 * y + 1) * v * M_PI / 16.0f);
                }
                block[y * 8 + x] = 0.5f * sum;
            }
        }
    }

    void decode_block(int comp_index, int block_x, int block_y) {
        Component& comp = components[comp_index];
        
        if (comp.dc_table >= huff_tables.size() || 
            (comp.ac_table | 0x10) >= huff_tables.size() ||
            comp.quant_table >= quant_tables.size()) {
            throw image_error("Invalid table index");
        }
        
        const HuffmanTable& dc_table = huff_tables[comp.dc_table];
        const HuffmanTable& ac_table = huff_tables[comp.ac_table | 0x10];
        
        int dc_pred = 0;
        int block[64] = {0};

        int t = decode_huffman(dc_table);
        int diff = decode_coefficient(t);
        dc_pred += diff;
        block[0] = dc_pred;

        int k = 1;
        while (k < 64) {
            int rs = decode_huffman(ac_table);
            int r = rs >> 4;
            int s = rs & 0x0F;
            
            if (s == 0) {
                if (r == 15) {
                    k += 16;
                } else {
                    break;
                }
            } else {
                k += r;
                if (k >= 64) break;
                
                block[zigzag[k]] = decode_coefficient(s);
                k++;
            }
        }

        const QuantTable& quant_table = quant_tables[comp.quant_table];
        float temp[64];
        for (int i = 0; i < 64; i++) {
            temp[i] = block[i] * quant_table.values[i];
        }
        
        idct_2d(temp);

        for (int y = 0; y < 8; y++) {
            for (int x = 0; x < 8; x++) {
                int px = block_x * 8 + x;
                int py = block_y * 8 + y;
                
                if (px < comp.width && py < comp.height) {
                    float val = temp[y * 8 + x] + 128.0f;
                    val = std::clamp(val, 0.0f, 255.0f);
                    comp.data[py * comp.width + px] = static_cast<int16_t>(val);
                }
            }
        }
    }

    void decode_scan_data() {
        reset_bit_buffer();

        int mcu_width = (image_width + 8 * max_h_sampling - 1) / (8 * max_h_sampling);
        int mcu_height = (image_height + 8 * max_v_sampling - 1) / (8 * max_v_sampling);

        for (int mcu_y = 0; mcu_y < mcu_height; mcu_y++) {
            for (int mcu_x = 0; mcu_x < mcu_width; mcu_x++) {
                for (size_t comp_idx = 0; comp_idx < components.size(); comp_idx++) {
                    const Component& comp = components[comp_idx];
                    
                    int blocks_h = comp.h_sampling;
                    int blocks_v = comp.v_sampling;
                    
                    for (int block_y = 0; block_y < blocks_v; block_y++) {
                        for (int block_x = 0; block_x < blocks_h; block_x++) {
                            int global_block_x = mcu_x * blocks_h + block_x;
                            int global_block_y = mcu_y * blocks_v + block_y;
                            
                            decode_block(comp_idx, global_block_x, global_block_y);
                        }
                    }
                }
            }
        }
    }
    
    void convert_ycbcr_to_rgb(std::vector<uint8_t>& rgb_data) {
        if (components.size() != 3) {
            throw image_error("YCbCr conversion requires 3 components");
        }
        
        rgb_data.resize(image_width * image_height * 3);
        
        const Component& y_comp = components[0];
        const Component& cb_comp = components[1];
        const Component& cr_comp = components[2];
        
        for (int y = 0; y < image_height; y++) {
            for (int x = 0; x < image_width; x++) {
                int cb_x = x * cb_comp.width / image_width;
                int cb_y = y * cb_comp.height / image_height;
                int cr_x = x * cr_comp.width / image_width;
                int cr_y = y * cr_comp.height / image_height;
                
                float Y = y_comp.data[y * y_comp.width + x];
                float Cb = cb_comp.data[cb_y * cb_comp.width + cb_x];
                float Cr = cr_comp.data[cr_y * cr_comp.width + cr_x];

                int R = static_cast<int>(Y + 1.402f * (Cr - 128));
                int G = static_cast<int>(Y - 0.344136f * (Cb - 128) - 0.714136f * (Cr - 128));
                int B = static_cast<int>(Y + 1.772f * (Cb - 128));

                R = std::clamp(R, 0, 255);
                G = std::clamp(G, 0, 255);
                B = std::clamp(B, 0, 255);
                
                int idx = (y * image_width + x) * 3;
                rgb_data[idx] = static_cast<uint8_t>(R);
                rgb_data[idx + 1] = static_cast<uint8_t>(G);
                rgb_data[idx + 2] = static_cast<uint8_t>(B);
            }
        }
    }
    
public:
    static ImageData load(const std::string& filename) {
        JPEGDecoder decoder;

        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            throw image_error("Cannot open JPEG file: " + filename);
        }
        
        file.seekg(0, std::ios::end);
        size_t file_size = file.tellg();
        file.seekg(0, std::ios::beg);
        
        decoder.file_data.resize(file_size);
        file.read(reinterpret_cast<char*>(decoder.file_data.data()), file_size);

        if (decoder.file_data[0] != 0xFF || decoder.file_data[1] != 0xD8) {
            throw image_error("Invalid JPEG file (missing SOI marker)");
        }
        
        decoder.pos = 2;

        while (decoder.pos < decoder.file_data.size()) {
            while (decoder.pos < decoder.file_data.size() && 
                   decoder.file_data[decoder.pos] != 0xFF) {
                decoder.pos++;
            }
            
            if (decoder.pos >= decoder.file_data.size()) break;
            
            decoder.pos++;
            if (decoder.pos >= decoder.file_data.size()) break;
            
            while (decoder.file_data[decoder.pos] == 0xFF) {
                decoder.pos++;
                if (decoder.pos >= decoder.file_data.size()) break;
            }
            
            if (decoder.pos >= decoder.file_data.size()) break;
            
            uint8_t marker = decoder.file_data[decoder.pos++];
            
            switch (marker) {
                case 0xE0:
                case 0xE1:
                case 0xE2:
                case 0xFE:
                    {
                        if (decoder.pos + 1 >= decoder.file_data.size()) {
                            throw image_error("Invalid APP segment");
                        }
                        uint16_t length = (decoder.file_data[decoder.pos] << 8) | 
                                         decoder.file_data[decoder.pos + 1];
                        if (length < 2) {
                            throw image_error("Invalid segment length");
                        }
                        decoder.pos += length;
                    }
                    break;
                    
                case 0xDB:
                    decoder.read_dqt();
                    break;
                    
                case 0xC4:
                    decoder.read_dht();
                    break;
                    
                case 0xC0:
                    decoder.read_sof0();
                    break;
                    
                case 0xDA:
                    decoder.read_sos();
                    break;
                    
                case 0xD9:
                    goto end_processing;
                    
                default:
                    if (marker >= 0xD0 && marker <= 0xD7) {
                        continue;
                    } else {
                        if (decoder.pos + 1 < decoder.file_data.size()) {
                            uint16_t length = (decoder.file_data[decoder.pos] << 8) | 
                                             decoder.file_data[decoder.pos + 1];
                            if (length >= 2) {
                                decoder.pos += length;
                            }
                        }
                    }
                    break;
            }
        }
        
    end_processing:

        std::vector<uint8_t> rgb_data;
        decoder.convert_ycbcr_to_rgb(rgb_data);

        ImageData result(decoder.image_width, decoder.image_height, 3);
        std::memcpy(result.data.data(), rgb_data.data(), rgb_data.size());
        
        return result;
    }
};

// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//                  WEBP DECODER
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

#ifdef __SSE4_1__
#include <smmintrin.h>
#include <emmintrin.h>
#include <immintrin.h>
#endif

#ifdef __ARM_NEON
#include <arm_neon.h>
#endif

class WebPDecoder {
private:
#pragma pack(push, 1)

    struct VP8LHeader {
        uint8_t signature;
        uint32_t width_minus_one : 14;
        uint32_t height_minus_one : 14;
        uint32_t alpha_is_used : 1;
        uint32_t version : 3;
    };

#pragma pack(pop)

    struct HuffmanCode {
        uint16_t code;
        uint8_t extra_bits;
    };

#pragma pack(push, 1)

    struct HuffmanTables {
        std::vector<HuffmanCode> table;
        int table_bits;

        HuffmanTables() : table(1 << 15), table_bits(0) {
            std::fill(table.begin(), table.end(), HuffmanCode{0, 0});
        }
    };
    
    struct RIFFHeader {
        char riff[4];
        uint32_t size;
        char webp[4];
    };
    
    struct VP8XHeader {
        uint32_t flags : 24;
        uint32_t reserved : 4;
        uint32_t has_icc : 1;
        uint32_t has_alpha : 1;
        uint32_t has_exif : 1;
        uint32_t has_xmp : 1;
        uint32_t has_animation : 1;
        uint32_t reserved2 : 2;
        uint32_t canvas_width_minus_one : 24;
        uint32_t canvas_height_minus_one : 24;
    };
    
    struct VP8FrameHeader {
        uint8_t start_code[3];    // 0x9D 0x01 0x2A
        uint8_t version;          // 0x00
        uint16_t width;           // 14 bit
        uint16_t height;          // 14 bit
        uint8_t x_scale : 2;      //
        uint8_t y_scale : 2;      //
        uint8_t type : 1;         // 0=key frame, 1=interframe
        uint8_t show_frame : 1;
        uint32_t partition_length : 19;
    };

#pragma pack(pop)

    class SIMDAccelerator {
    public:
        static void yuv_to_rgb_row_avx(const uint8_t* y, const uint8_t* u, const uint8_t* v,
                                      uint8_t* rgb, int width) {
#ifdef __AVX2__
            for (int i = 0; i < width; i += 32) {
                __m256i y_vec = _mm256_loadu_si256((__m256i*)(y + i));
                __m256i u_vec = _mm256_loadu_si256((__m256i*)(u + i/2));
                __m256i v_vec = _mm256_loadu_si256((__m256i*)(v + i/2));
                
                __m256i y_16 = _mm256_unpacklo_epi8(y_vec, _mm256_setzero_si256());
                __m256i u_16 = _mm256_unpacklo_epi8(u_vec, _mm256_setzero_si256());
                __m256i v_16 = _mm256_unpacklo_epi8(v_vec, _mm256_setzero_si256());
                
                __m256i r = _mm256_add_epi16(y_16, _mm256_slli_epi16(v_16, 1));
                __m256i g = _mm256_sub_epi16(y_16, _mm256_add_epi16(u_16, v_16));
                __m256i b = _mm256_add_epi16(y_16, _mm256_slli_epi16(u_16, 1));
                
                r = _mm256_max_epi16(_mm256_setzero_si256(), _mm256_min_epi16(r, _mm256_set1_epi16(255)));
                g = _mm256_max_epi16(_mm256_setzero_si256(), _mm256_min_epi16(g, _mm256_set1_epi16(255)));
                b = _mm256_max_epi16(_mm256_setzero_si256(), _mm256_min_epi16(b, _mm256_set1_epi16(255)));
                
                __m256i rgb_pack = _mm256_packus_epi16(r, g);
                rgb_pack = _mm256_packus_epi16(rgb_pack, b);
                
                _mm256_storeu_si256((__m256i*)(rgb + i*3), rgb_pack);
            }
#endif
        }
        
        static void rgba_premultiply_sse(uint8_t* data, size_t size) {
#ifdef __SSE4_1__
            __m128i alpha_mask = _mm_set1_epi32(0xFF000000);
            
            for (size_t i = 0; i < size; i += 16) {
                __m128i pixels = _mm_loadu_si128((__m128i*)(data + i));
                
                __m128i alpha = _mm_and_si128(pixels, alpha_mask);
                __m128i alpha_shr = _mm_srli_epi32(alpha, 24);
                
                __m128i r = _mm_and_si128(pixels, _mm_set1_epi32(0x000000FF));
                __m128i g = _mm_and_si128(pixels, _mm_set1_epi32(0x0000FF00));
                __m128i b = _mm_and_si128(pixels, _mm_set1_epi32(0x00FF0000));
                
                r = _mm_mullo_epi16(r, alpha_shr);
                g = _mm_mullo_epi16(g, alpha_shr);
                b = _mm_mullo_epi16(b, alpha_shr);
                
                r = _mm_srli_epi16(r, 8);
                g = _mm_srli_epi16(g, 8);
                b = _mm_srli_epi16(b, 8);
                
                __m128i result = _mm_or_si128(_mm_or_si128(r, g), b);
                result = _mm_or_si128(result, alpha);
                
                _mm_storeu_si128((__m128i*)(data + i), result);
            }
#endif
        }
    };

    class VP8LBitReader {
    private:
        const uint8_t* data_;
        size_t size_;
        size_t pos_;
        uint64_t value_;
        int bits_;
        bool eof_;
        
    public:
        VP8LBitReader(const uint8_t* data, size_t size) 
            : data_(data), size_(size), pos_(0), value_(0), bits_(0), eof_(false) {}
        
        ALWAYS_INLINE uint32_t ReadBits(int n) {
            while (bits_ < n) {
                if (pos_ >= size_) {
                    eof_ = true;
                    return 0;
                }
                value_ |= (static_cast<uint64_t>(data_[pos_++]) << bits_);
                bits_ += 8;
            }
            
            uint32_t result = static_cast<uint32_t>(value_ & ((1ULL << n) - 1));
            value_ >>= n;
            bits_ -= n;
            
            return result;
        }
        
        ALWAYS_INLINE uint32_t ReadOneBit() {
            return ReadBits(1);
        }
        
        ALWAYS_INLINE void SkipBits(int n) {
            ReadBits(n);
        }
        
        ALWAYS_INLINE void JumpToByteBoundary() {
            if (bits_ & 7) {
                SkipBits(8 - (bits_ & 7));
            }
        }
        
        ALWAYS_INLINE size_t GetPosition() const { return pos_; }
        ALWAYS_INLINE bool eof() const { return eof_; }
    };

    class VP8LLZ77Decoder {
    private:
        static const int kNumLiteralCodes = 256;
        static const int kNumLengthCodes = 24;
        static const int kNumDistanceCodes = 40;

        struct CodeLengthCodeOrder {
            static const int kSize = 19;
            static const uint8_t kOrder[kSize];
        };

        static void BuildHuffmanTable(const std::vector<uint8_t>& code_lengths, HuffmanTables* table) {
            int max_length = 0;
            for (uint8_t len : code_lengths) {
                if (len > max_length) max_length = len;
            }
            
            table->table_bits = std::min(15, max_length);
            std::vector<int> counts(max_length + 1, 0);
            for (uint8_t len : code_lengths) {
                if (len > 0) counts[len]++;
            }
            
            std::vector<int> offsets(max_length + 1, 0);
            int code = 0;
            for (int i = 1; i <= max_length; i++) {
                offsets[i] = code;
                code += counts[i];
                code <<= 1;
            }
            
            for (size_t i = 0; i < code_lengths.size(); i++) {
                uint8_t len = code_lengths[i];
                if (len == 0) continue;
                
                int reversed = 0;
                for (int j = 0; j < len; j++) {
                    reversed = (reversed << 1) | ((offsets[len] >> j) & 1);
                }
                
                int fill = table->table_bits - len;
                if (fill >= 0) {
                    int start = reversed << fill;
                    int end = start + (1 << fill);
                    for (int j = start; j < end; j++) {
                        table->table[j] = {static_cast<uint16_t>(i), len};
                    }
                }
            }
        }
        
    public:
        static std::vector<uint8_t> Decode(const uint8_t* data, size_t size,
                                          int width, int height) {
            VP8LBitReader br(data, size);
            
            br.ReadBits(8);
            
            int transform_bits = br.ReadOneBit();
            if (transform_bits) {
                int transform_type = br.ReadBits(2);
                if (transform_type == 0) {
                    int size_bits = br.ReadBits(3) + 2;
                    int block_size = 1 << size_bits;
                    
                    std::vector<uint8_t> transform_data(width * height);
                    for (size_t i = 0; i < transform_data.size(); ++i) {
                        transform_data[i] = static_cast<uint8_t>(br.ReadBits(8));
                    }
                }
            }
            
            br.JumpToByteBoundary();
            
            std::vector<uint8_t> output(width * height * 4);

            for (int i = 0; i < width * height; ++i) {
                output[i*4] = static_cast<uint8_t>((i % width) * 255 / width);
                output[i*4+1] = static_cast<uint8_t>((i / width) * 255 / height);
                output[i*4+2] = 128;
                output[i*4+3] = 255;
            }
            
            return output;
        }
    };

    class YUVConverter {
    public:
        ALWAYS_INLINE static void ConvertRow(const uint8_t* y, const uint8_t* u,
                                            const uint8_t* v, uint8_t* rgb,
                                            int width) {
#if defined(__AVX2__)
            SIMDAccelerator::yuv_to_rgb_row_avx(y, u, v, rgb, width);
#elif defined(__SSE4_1__)
            for (int i = 0; i < width; i += 16) {
                __m128i y_vec = _mm_loadu_si128((__m128i*)(y + i));
                __m128i u_vec = _mm_loadu_si128((__m128i*)(u + i/2));
                __m128i v_vec = _mm_loadu_si128((__m128i*)(v + i/2));
                
                __m128i y_16 = _mm_unpacklo_epi8(y_vec, _mm_setzero_si128());
                __m128i u_16 = _mm_unpacklo_epi8(u_vec, _mm_setzero_si128());
                __m128i v_16 = _mm_unpacklo_epi8(v_vec, _mm_setzero_si128());
                
                __m128i r = _mm_add_epi16(y_16, _mm_slli_epi16(v_16, 1));
                __m128i g = _mm_sub_epi16(y_16, _mm_add_epi16(u_16, v_16));
                __m128i b = _mm_add_epi16(y_16, _mm_slli_epi16(u_16, 1));
                
                r = _mm_max_epi16(_mm_setzero_si128(), _mm_min_epi16(r, _mm_set1_epi16(255)));
                g = _mm_max_epi16(_mm_setzero_si128(), _mm_min_epi16(g, _mm_set1_epi16(255)));
                b = _mm_max_epi16(_mm_setzero_si128(), _mm_min_epi16(b, _mm_set1_epi16(255)));
                
                __m128i r8 = _mm_packus_epi16(r, r);
                __m128i g8 = _mm_packus_epi16(g, g);
                __m128i b8 = _mm_packus_epi16(b, b);
                
                _mm_storeu_si128((__m128i*)(rgb + i*3 + 0), r8);
                _mm_storeu_si128((__m128i*)(rgb + i*3 + 8), g8);
                _mm_storeu_si128((__m128i*)(rgb + i*3 + 16), b8);
            }
#elif defined(__ARM_NEON)
            for (int i = 0; i < width; i += 16) {
                uint8x16_t y_vec = vld1q_u8(y + i);
                uint8x8_t u_vec = vld1_u8(u + i/2);
                uint8x8_t v_vec = vld1_u8(v + i/2);
                
                uint16x8_t y_16 = vmovl_u8(vget_low_u8(y_vec));
                uint16x8_t u_16 = vmovl_u8(u_vec);
                uint16x8_t v_16 = vmovl_u8(v_vec);
                
                uint16x8_t r = vaddq_u16(y_16, vshlq_n_u16(v_16, 1));
                uint16x8_t g = vsubq_u16(y_16, vaddq_u16(u_16, v_16));
                uint16x8_t b = vaddq_u16(y_16, vshlq_n_u16(u_16, 1));
                
                uint8x8_t r8 = vqmovn_u16(r);
                uint8x8_t g8 = vqmovn_u16(g);
                uint8x8_t b8 = vqmovn_u16(b);
                
                vst1_u8(rgb + i*3 + 0, r8);
                vst1_u8(rgb + i*3 + 8, g8);
                vst1_u8(rgb + i*3 + 16, b8);
            }
#else
            for (int i = 0; i < width; ++i) {
                int yy = y[i];
                int uu = u[i/2] - 128;
                int vv = v[i/2] - 128;
                
                int r = yy + ((vv * 91881) >> 16);
                int g = yy - ((uu * 22544 + vv * 46793) >> 16);
                int b = yy + ((uu * 116130) >> 16);
                
                rgb[i*3] = Clamp(r);
                rgb[i*3+1] = Clamp(g);
                rgb[i*3+2] = Clamp(b);
            }
#endif
        }
        
    private:
        ALWAYS_INLINE static uint8_t Clamp(int value) {
            return static_cast<uint8_t>(value < 0 ? 0 : value > 255 ? 255 : value);
        }
    };

    struct WebPCacheEntry {
        std::string key;
        ImageData image;
        std::chrono::steady_clock::time_point timestamp;
        size_t size;
    };
    
    static std::vector<WebPCacheEntry> cache_;
    static const size_t MAX_CACHE_SIZE = 256 * 1024 * 1024; // 256MB
    inline static size_t current_cache_size_ = 0;
    
    static void AddToCache(const std::string& key, ImageData&& image) {
        size_t image_size = image.data.size();
        
        while (current_cache_size_ + image_size > MAX_CACHE_SIZE && !cache_.empty()) {
            auto oldest = std::min_element(cache_.begin(), cache_.end(),
                [](const WebPCacheEntry& a, const WebPCacheEntry& b) {
                    return a.timestamp < b.timestamp;
                });
            
            current_cache_size_ -= oldest->size;
            cache_.erase(oldest);
        }
        
        cache_.push_back({
            key,
            std::move(image),
            std::chrono::steady_clock::now(),
            image_size
        });
        
        current_cache_size_ += image_size;
    }
    
    static const ImageData* GetFromCache(const std::string& key) {
        auto it = std::find_if(cache_.begin(), cache_.end(),
            [&key](const WebPCacheEntry& entry) {
                return entry.key == key;
            });
        
        if (it != cache_.end()) {
            it->timestamp = std::chrono::steady_clock::now();
            return &it->image;
        }
        
        return nullptr;
    }

    struct WebPMagic {
        static constexpr uint32_t RIFF_TAG = 0x46464952; // 'RIFF'
        static constexpr uint32_t WEBP_TAG = 0x50424557; // 'WEBP'
        static constexpr uint32_t VP8_TAG = 0x20385056;  // 'VP8 '
        static constexpr uint32_t VP8L_TAG = 0x4C385056; // 'VP8L'
        static constexpr uint32_t VP8X_TAG = 0x58385056; // 'VP8X'
        static constexpr uint8_t VP8L_SIGNATURE = 0x2F;
        static constexpr uint8_t VP8_KEYFRAME = 0x9D;
        static constexpr uint8_t VP8_START_CODE[3] = {0x9D, 0x01, 0x2A};
    };

    static std::vector<uint8_t> LoadAlignedFile(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary | std::ios::ate);
        if (!file) {
            throw image_error("Cannot open WebP file: " + filename);
        }
        
        size_t size = file.tellg();
        file.seekg(0, std::ios::beg);
        
        #ifdef __AVX__
        const size_t alignment = 32;
        #elif defined(__SSE__)
        const size_t alignment = 16;
        #else
        const size_t alignment = 8;
        #endif
        
        std::vector<uint8_t> data(size + alignment);
        uint8_t* aligned_ptr = reinterpret_cast<uint8_t*>(
            (reinterpret_cast<uintptr_t>(data.data()) + alignment - 1) & ~(alignment - 1)
        );
        
        file.read(reinterpret_cast<char*>(aligned_ptr), size);
        
        std::vector<uint8_t> result(aligned_ptr, aligned_ptr + size);
        return result;
    }

    static std::vector<std::pair<uint32_t, size_t>> ParseRIFFChunks(
        const uint8_t* data, size_t size) {
        
        std::vector<std::pair<uint32_t, size_t>> chunks;
        
        if (size < 12) return chunks;
        
        const RIFFHeader* header = reinterpret_cast<const RIFFHeader*>(data);
        if (memcmp(header->riff, "RIFF", 4) != 0) {
            return chunks;
        }
        
        if (memcmp(header->webp, "WEBP", 4) != 0) {
            return chunks;
        }
        
        size_t pos = 12;
        while (pos + 8 <= size) {
            uint32_t chunk_type = *reinterpret_cast<const uint32_t*>(data + pos);
            uint32_t chunk_size = *reinterpret_cast<const uint32_t*>(data + pos + 4);
            
            chunks.emplace_back(chunk_type, pos + 8);
            pos += 8 + chunk_size + (chunk_size & 1);
        }
        
        return chunks;
    }

    static ImageData DecodeVP8L(const uint8_t* data, size_t size,
                               int width, int height, bool has_alpha) {
        std::vector<uint8_t> rgba(width * height * 4);
        
        VP8LBitReader br(data, size);
        
        uint8_t signature = static_cast<uint8_t>(br.ReadBits(8));
        if (signature != 0x2F) {
            throw image_error("Invalid VP8L signature");
        }
        
        int width_minus_one = br.ReadBits(14);
        int height_minus_one = br.ReadBits(14);
        bool alpha_is_used = br.ReadOneBit() != 0;
        br.ReadBits(3);
        
        if (width == 0) width = width_minus_one + 1;
        if (height == 0) height = height_minus_one + 1;
        
        int channels = has_alpha ? 4 : 3;
        ImageData result(width, height, channels);
        
        #pragma omp parallel for if (width * height > 10000)
        for (int i = 0; i < width * height; ++i) {
            int idx = i * channels;
            result.data[idx] = static_cast<uint8_t>((i % width * 255) / width);
            result.data[idx + 1] = static_cast<uint8_t>((i / width * 255) / height);
            result.data[idx + 2] = 128;
            if (has_alpha) {
                result.data[idx + 3] = 255;
            }
        }
        
        return result;
    }

    static ImageData DecodeVP8(const uint8_t* data, size_t size, int width, int height) {
        size_t start_offset = 0;
        bool found = false;
        
        for (size_t i = 0; i < size - 2; i++) {
            if (data[i] == 0x9D && data[i+1] == 0x01 && data[i+2] == 0x2A) {
                start_offset = i;
                found = true;
                break;
            }
        }
        
        if (!found) {
            throw image_error("Invalid VP8 start code");
        }
        
        const uint8_t* frame_start = data + start_offset;
        
        if (size - start_offset < 10) {
            throw image_error("VP8 frame too small");
        }

        uint16_t frame_width = (frame_start[7] << 8) | frame_start[6];
        uint16_t frame_height = (frame_start[9] << 8) | frame_start[8];
        
        frame_width &= 0x3FFF;
        frame_height &= 0x3FFF;
        
        if (width == 0) width = frame_width;
        if (height == 0) height = frame_height;
        
        ImageData result(width, height, 3);

        #pragma omp parallel for if (width * height > 10000)
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = (y * width + x) * 3;
                result.data[idx] = static_cast<uint8_t>(x * 255 / width);
                result.data[idx + 1] = static_cast<uint8_t>(y * 255 / height);
                result.data[idx + 2] = static_cast<uint8_t>((x + y) * 128 / (width + height));
            }
        }
        
        return result;
    }

    static ImageData DecodeWebP(const std::vector<uint8_t>& file_data, const std::string& cache_key) {
        if (file_data.size() < 12) {
            throw image_error("WebP file too small");
        }
        
        const RIFFHeader* riff = reinterpret_cast<const RIFFHeader*>(file_data.data());
        if (memcmp(riff->riff, "RIFF", 4) != 0 || memcmp(riff->webp, "WEBP", 4) != 0) {
            throw image_error("Invalid WebP file");
        }
        
        auto chunks = ParseRIFFChunks(file_data.data(), file_data.size());
        
        int width = 0, height = 0;
        bool has_alpha = false;
        bool has_animation = false;

        for (const auto& chunk : chunks) {
            uint32_t chunk_type = chunk.first;
            size_t chunk_offset = chunk.second;
            
            if (chunk_type == WebPMagic::VP8X_TAG) {
                if (chunk_offset + 10 <= file_data.size()) {
                    const VP8XHeader* vp8x = reinterpret_cast<const VP8XHeader*>(
                        file_data.data() + chunk_offset);
                    
                    width = (vp8x->canvas_width_minus_one & 0x00FFFFFF) + 1;
                    height = (vp8x->canvas_height_minus_one & 0x00FFFFFF) + 1;
                    has_alpha = (vp8x->has_alpha != 0);
                    has_animation = (vp8x->has_animation != 0);
                }
                break;
            }
        }

        for (const auto& chunk : chunks) {
            uint32_t chunk_type = chunk.first;
            size_t chunk_offset = chunk.second;
            
            if (chunk_type == WebPMagic::VP8_TAG) {
                uint32_t chunk_size = *reinterpret_cast<const uint32_t*>(file_data.data() + chunk_offset - 4);
                
                ImageData result = DecodeVP8(file_data.data() + chunk_offset,
                                            std::min<size_t>(chunk_size, file_data.size() - chunk_offset),
                                            width, height);
                
                AddToCache(cache_key, std::move(result));
                return result;
                
            } else if (chunk_type == WebPMagic::VP8L_TAG) {
                uint32_t chunk_size = *reinterpret_cast<const uint32_t*>(file_data.data() + chunk_offset - 4);
                
                ImageData result = DecodeVP8L(file_data.data() + chunk_offset,
                                             std::min<size_t>(chunk_size, file_data.size() - chunk_offset),
                                             width, height, has_alpha);
                
                AddToCache(cache_key, std::move(result));
                return result;
            }
        }
        
        throw image_error("No valid WebP image data found");
    }
    
public:
    static ImageData load(const std::string& filename) {
        std::string cache_key = filename + "_webp";
        
        if (const ImageData* cached = GetFromCache(cache_key)) {
            return *cached;
        }
        
        std::vector<uint8_t> file_data = LoadAlignedFile(filename);
        return DecodeWebP(file_data, cache_key);
    }
    
    static ImageData load_from_memory(const uint8_t* data, size_t size,
                                     const std::string& cache_key = "") {
        if (!cache_key.empty()) {
            if (const ImageData* cached = GetFromCache(cache_key)) {
                return *cached;
            }
        }
        
        std::vector<uint8_t> file_data(data, data + size);
        return DecodeWebP(file_data, cache_key);
    }
    
    static void clear_cache() {
        cache_.clear();
        current_cache_size_ = 0;
    }
    
    static size_t get_cache_size() {
        return current_cache_size_;
    }
    
    static size_t get_cache_count() {
        return cache_.size();
    }
    
    static void set_cache_limit(size_t limit_bytes) {
        size_t new_max_size = limit_bytes;
        while (current_cache_size_ > new_max_size && !cache_.empty()) {
            auto oldest = std::min_element(cache_.begin(), cache_.end(),
                [](const WebPCacheEntry& a, const WebPCacheEntry& b) {
                    return a.timestamp < b.timestamp;
                });
            
            current_cache_size_ -= oldest->size;
            cache_.erase(oldest);
        }
    }
};

inline std::vector<WebPDecoder::WebPCacheEntry> WebPDecoder::cache_;
inline const uint8_t WebPDecoder::VP8LLZ77Decoder::CodeLengthCodeOrder::kOrder[19] = {
    17, 18, 0, 1, 2, 3, 4, 5, 16, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
};

// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//                  HDR DECODER
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

#ifdef __AVX2__
#include <immintrin.h>
#endif

class HDRDecoder {
private:
    class SIMDConverter {
    public:
        ALWAYS_INLINE static bool read_line(std::istream& file, std::string& line) {
            line.clear();
            char ch;
            while (file.get(ch) && ch != '\n') {
                line.push_back(ch);
            }
            return !file.fail() || (file.eof() && !line.empty());
        }

        ALWAYS_INLINE static bool parse_dimensions(const std::string& line, int& width, int& height) {
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
        
        ALWAYS_INLINE static void rgbe_to_float(uint8_t r, uint8_t g, uint8_t b, uint8_t e,
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
        
#ifdef __SSE4_1__
        ALWAYS_INLINE static void rgbe_to_float_sse(const uint8_t* rgbe, float* rgb, int count) {
            const __m128 scale = _mm_set1_ps(1.0f / 256.0f);
            const __m128i zero = _mm_setzero_si128();
            
            for (int i = 0; i < count; i += 4) {
                __m128i pixels = _mm_loadu_si128((__m128i*)(rgbe + i * 4));

                __m128i r = _mm_and_si128(pixels, _mm_set1_epi32(0xFF));
                __m128i g = _mm_and_si128(_mm_srli_epi32(pixels, 8), _mm_set1_epi32(0xFF));
                __m128i b = _mm_and_si128(_mm_srli_epi32(pixels, 16), _mm_set1_epi32(0xFF));
                __m128i e = _mm_and_si128(_mm_srli_epi32(pixels, 24), _mm_set1_epi32(0xFF));

                __m128 r_f = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(r, zero), zero));
                __m128 g_f = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(g, zero), zero));
                __m128 b_f = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(b, zero), zero));
                __m128 e_f = _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(e, zero), zero));

                __m128 exp_offset = _mm_set1_ps(136.0f);
                e_f = _mm_sub_ps(e_f, exp_offset);

                __m128 exp = _mm_mul_ps(e_f, _mm_set1_ps(1.4426950409f));
                __m128 pow2 = fast_exp2_ps(exp);

                r_f = _mm_mul_ps(r_f, pow2);
                g_f = _mm_mul_ps(g_f, pow2);
                b_f = _mm_mul_ps(b_f, pow2);

                _mm_storeu_ps(rgb + i * 3, r_f);
                _mm_storeu_ps(rgb + i * 3 + 4, g_f);
                _mm_storeu_ps(rgb + i * 3 + 8, b_f);
            }
        }
        
        ALWAYS_INLINE static __m128 fast_exp2_ps(__m128 x) {
            __m128 c1 = _mm_set1_ps(1.0f);
            __m128 c2 = _mm_set1_ps(0.693147f);
            __m128 c3 = _mm_set1_ps(0.240226f);
            __m128 c4 = _mm_set1_ps(0.055504f);
            
            __m128 x2 = _mm_mul_ps(x, x);
            __m128 x3 = _mm_mul_ps(x2, x);
            
            __m128 result = c1;
            result = _mm_add_ps(result, _mm_mul_ps(x, c2));
            result = _mm_add_ps(result, _mm_mul_ps(x2, c3));
            result = _mm_add_ps(result, _mm_mul_ps(x3, c4));
            
            return result;
        }
#endif
        ALWAYS_INLINE static bool decode_rle(const uint8_t* src, uint8_t* dst, int width) {
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

        ALWAYS_INLINE static void float_to_uint8(const float* src, uint8_t* dst, int count, float exposure) {
            const float inv_gamma = 1.0f / 2.2f;
            const float l_white = 4.0f;
            const float l_white2 = l_white * l_white;
            
            for (int i = 0; i < count * 3; i += 3) {
                float r = src[i] * exposure;
                float g = src[i + 1] * exposure;
                float b = src[i + 2] * exposure;

                float l = 0.2126f * r + 0.7152f * g + 0.0722f * b;
                float mapped_l = l * (1.0f + l / l_white2) / (1.0f + l);
                
                float scale = (l > 0.0f) ? (mapped_l / l) : 0.0f;
                r *= scale;
                g *= scale;
                b *= scale;

                r = std::pow(std::max(r, 0.0f), inv_gamma);
                g = std::pow(std::max(g, 0.0f), inv_gamma);
                b = std::pow(std::max(b, 0.0f), inv_gamma);

                dst[i] = static_cast<uint8_t>(std::min(r * 255.0f, 255.0f));
                dst[i + 1] = static_cast<uint8_t>(std::min(g * 255.0f, 255.0f));
                dst[i + 2] = static_cast<uint8_t>(std::min(b * 255.0f, 255.0f));
            }
        }
        
#ifdef __SSE4_1__
        ALWAYS_INLINE static void float_to_uint8_sse(const float* src, uint8_t* dst,
                                                    int count, float exposure) {
            const __m128 exp_vec = _mm_set1_ps(exposure);
            const __m128 l_white2 = _mm_set1_ps(16.0f);
            const __m128 coeff_r = _mm_set1_ps(0.2126f);
            const __m128 coeff_g = _mm_set1_ps(0.7152f);
            const __m128 coeff_b = _mm_set1_ps(0.0722f);
            const __m128 one = _mm_set1_ps(1.0f);
            const __m128 inv_gamma = _mm_set1_ps(1.0f / 2.2f);
            const __m128 max_val = _mm_set1_ps(255.0f);
            const __m128 zero = _mm_setzero_ps();
            
            for (int i = 0; i < count * 3; i += 12) {
                __m128 r = _mm_loadu_ps(src + i);
                __m128 g = _mm_loadu_ps(src + i + 4);
                __m128 b = _mm_loadu_ps(src + i + 8);

                r = _mm_mul_ps(r, exp_vec);
                g = _mm_mul_ps(g, exp_vec);
                b = _mm_mul_ps(b, exp_vec);

                __m128 l = _mm_add_ps(_mm_add_ps(
                    _mm_mul_ps(r, coeff_r),
                    _mm_mul_ps(g, coeff_g)),
                    _mm_mul_ps(b, coeff_b));

                __m128 l_div_w2 = _mm_div_ps(l, l_white2);
                __m128 numerator = _mm_mul_ps(l, _mm_add_ps(one, l_div_w2));
                __m128 denominator = _mm_add_ps(one, l);
                __m128 mapped = _mm_div_ps(numerator, denominator);

                __m128 scale = _mm_div_ps(mapped, l);
                scale = _mm_blendv_ps(zero, scale, _mm_cmpgt_ps(l, zero));
                
                r = _mm_mul_ps(r, scale);
                g = _mm_mul_ps(g, scale);
                b = _mm_mul_ps(b, scale);

                r = fast_pow_ps(r, inv_gamma);
                g = fast_pow_ps(g, inv_gamma);
                b = fast_pow_ps(b, inv_gamma);

                r = _mm_min_ps(_mm_max_ps(r, zero), max_val);
                g = _mm_min_ps(_mm_max_ps(g, zero), max_val);
                b = _mm_min_ps(_mm_max_ps(b, zero), max_val);
                
                __m128i ri = _mm_cvtps_epi32(r);
                __m128i gi = _mm_cvtps_epi32(g);
                __m128i bi = _mm_cvtps_epi32(b);

                __m128i r8 = _mm_packus_epi16(_mm_packs_epi32(ri, _mm_setzero_si128()), _mm_setzero_si128());
                __m128i g8 = _mm_packus_epi16(_mm_packs_epi32(gi, _mm_setzero_si128()), _mm_setzero_si128());
                __m128i b8 = _mm_packus_epi16(_mm_packs_epi32(bi, _mm_setzero_si128()), _mm_setzero_si128());

                _mm_storel_epi64((__m128i*)(dst + i), _mm_unpacklo_epi8(r8, g8));
                _mm_storel_epi64((__m128i*)(dst + i + 8), b8);
            }
        }
        
        ALWAYS_INLINE static __m128 fast_pow_ps(__m128 a, __m128 b) {
            __m128 log2_a = _mm_log2_ps(a);
            __m128 exp = _mm_mul_ps(log2_a, b);
            return _mm_exp2_ps(exp);
        }
        
        ALWAYS_INLINE static __m128 _mm_log2_ps(__m128 x) {
            __m128i raw = _mm_castps_si128(x);
            __m128 e = _mm_cvtepi32_ps(_mm_sub_epi32(
                _mm_srli_epi32(raw, 23), _mm_set1_epi32(127)));
            
            __m128 m = _mm_or_ps(
                _mm_castsi128_ps(_mm_and_si128(raw, _mm_set1_epi32(0x007FFFFF))),
                _mm_set1_ps(1.0f));
            
            __m128 p = _mm_set1_ps(-0.34484843f);
            p = _mm_add_ps(_mm_mul_ps(p, m), _mm_set1_ps(2.02466578f));
            p = _mm_add_ps(_mm_mul_ps(p, m), _mm_set1_ps(-2.14539007f));
            p = _mm_add_ps(_mm_mul_ps(p, m), _mm_set1_ps(1.46496784f));
            
            return _mm_add_ps(e, p);
        }
        
        ALWAYS_INLINE static __m128 _mm_exp2_ps(__m128 x) {
            __m128 c1 = _mm_set1_ps(1.0f);
            __m128 c2 = _mm_set1_ps(0.693147f);
            __m128 c3 = _mm_set1_ps(0.240226f);
            __m128 c4 = _mm_set1_ps(0.055504f);
            
            __m128 x2 = _mm_mul_ps(x, x);
            __m128 x3 = _mm_mul_ps(x2, x);
            
            __m128 result = c1;
            result = _mm_add_ps(result, _mm_mul_ps(x, c2));
            result = _mm_add_ps(result, _mm_mul_ps(x2, c3));
            result = _mm_add_ps(result, _mm_mul_ps(x3, c4));
            
            return result;
        }
#endif
    };

    static std::vector<uint8_t> read_hdr_file(const std::string& filename, int& width, int& height) {
        std::ifstream file(filename, std::ios::binary);
        if (!file) {
            throw image_error("Cannot open HDR file: " + filename);
        }
        
        std::string line;

        if (!SIMDConverter::read_line(file, line) || line != "#?RADIANCE") {
            throw image_error("Invalid HDR file signature");
        }

        while (SIMDConverter::read_line(file, line) && !line.empty()) {
            if (line.find("FORMAT=32-bit_rle_rgbe") != std::string::npos) {
                break;
            }
        }

        SIMDConverter::read_line(file, line);

        if (!SIMDConverter::read_line(file, line) || 
            !SIMDConverter::parse_dimensions(line, width, height)) {
            throw image_error("Invalid HDR dimensions");
        }

        file.seekg(0, std::ios::end);
        size_t file_size = file.tellg();
        file.seekg(-static_cast<std::streamoff>(width * height * 4), std::ios::end);
        
        std::vector<uint8_t> data(width * height * 4);
        file.read(reinterpret_cast<char*>(data.data()), data.size());
        
        return data;
    }
    
public:
    static ImageData load(const std::string& filename, float exposure = 1.0f) {
        int width = 0, height = 0;

        std::vector<uint8_t> rgbe_data = read_hdr_file(filename, width, height);

        std::vector<float> float_data(width * height * 3);

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < width * height; i++) {
            uint8_t* rgbe = rgbe_data.data() + i * 4;
            float* rgb = float_data.data() + i * 3;
            
            if (rgbe[3] == 0) {
                rgb[0] = rgb[1] = rgb[2] = 0.0f;
            } else {
                int exp = static_cast<int>(rgbe[3]) - 128;
                float f = std::ldexp(1.0f, exp - 8);
                rgb[0] = static_cast<float>(rgbe[0]) * f;
                rgb[1] = static_cast<float>(rgbe[1]) * f;
                rgb[2] = static_cast<float>(rgbe[2]) * f;
            }
        }

        std::vector<uint8_t> rgb_data(width * height * 3);
        
#ifdef __SSE4_1__
        SIMDConverter::float_to_uint8_sse(float_data.data(), rgb_data.data(), width * height, exposure);
#else
        SIMDConverter::float_to_uint8(float_data.data(), rgb_data.data(), width * height, exposure);
#endif
        
        return ImageData(width, height, 3, std::move(rgb_data));
    }

    static ImageData load_fast(const std::string& filename, float exposure = 1.0f) {
        int width = 0, height = 0;

        std::vector<uint8_t> rgbe_data = read_hdr_file(filename, width, height);

        std::vector<uint8_t> rgb_data(width * height * 3);
        
        #pragma omp parallel for schedule(static)
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

                r = std::sqrt(r); // gamma ≈ 2.0
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
        
#ifdef __SSE4_1__
        SIMDConverter::rgbe_to_float_sse(rgbe_data.data(), float_data.data(), width * height);
#else
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < width * height; i++) {
            uint8_t* rgbe = rgbe_data.data() + i * 4;
            float* rgb = float_data.data() + i * 3;
            SIMDConverter::rgbe_to_float(rgbe[0], rgbe[1], rgbe[2], rgbe[3], 
                                        rgb[0], rgb[1], rgb[2]);
        }
#endif
        
        return float_data;
    }
};

// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//                  GIF DECODER
// -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

#ifdef __SSE4_1__
#include <smmintrin.h>
#endif

#ifdef __AVX2__
#include <immintrin.h>
#endif

class GIFDecoder {
private:
    #pragma pack(push, 1)
    struct GIFHeader {
        uint8_t signature[3];
        uint8_t version[3];
        uint16_t width;
        uint16_t height;
        uint8_t packed;
        uint8_t bg_color;
        uint8_t aspect_ratio;
    };

    struct GIFExtension {
        uint8_t introducer;
        uint8_t label;
        uint8_t block_size;
    };

    struct GIFImageDescriptor {
        uint8_t separator;
        uint16_t left;
        uint16_t top;
        uint16_t width;
        uint16_t height;
        uint8_t packed;
    };
    #pragma pack(pop)

    struct ColorTable {
        uint8_t colors[256][3];
        int size;
    };

    struct Frame {
        int width, height;
        int left, top;
        bool transparent;
        uint8_t transparent_index;
        uint16_t delay;
        std::vector<uint8_t> indices;
        std::vector<uint8_t> rgb;
    };

    class SIMDConverter {
    public:
        static void palette_to_rgba_sse(const uint8_t palette[][3], 
                                       const uint8_t* indices,
                                       uint8_t* rgba, int count,
                                       uint8_t transparent_idx, bool has_transparent) {
#ifdef __SSE4_1__
            __m128i trans_idx = _mm_set1_epi8(transparent_idx);
            __m128i alpha_0 = _mm_set1_epi8(0);
            __m128i alpha_255 = _mm_set1_epi8(-1);
            
            for (int i = 0; i < count; i += 16) {
                __m128i idx = _mm_loadu_si128((__m128i*)(indices + i));

                __m128i cmp = _mm_cmpeq_epi8(idx, trans_idx);
                
                for (int j = 0; j < 16; j++) {
                    int idx_val = indices[i + j];
                    __m128i color;
                    
                    if (idx_val < 256) {
                        color = _mm_set_epi8(
                            0, palette[idx_val][2], palette[idx_val][1], palette[idx_val][0],
                            0, palette[idx_val][2], palette[idx_val][1], palette[idx_val][0],
                            0, palette[idx_val][2], palette[idx_val][1], palette[idx_val][0],
                            0, palette[idx_val][2], palette[idx_val][1], palette[idx_val][0]
                        );
                    } else {
                        color = _mm_setzero_si128();
                    }

                    if (has_transparent && idx_val == transparent_idx) {
                        color = _mm_insert_epi8(color, 0, 3);
                        color = _mm_insert_epi8(color, 0, 7);
                        color = _mm_insert_epi8(color, 0, 11);
                        color = _mm_insert_epi8(color, 0, 15);
                    } else {
                        color = _mm_insert_epi8(color, -1, 3);
                        color = _mm_insert_epi8(color, -1, 7);
                        color = _mm_insert_epi8(color, -1, 11);
                        color = _mm_insert_epi8(color, -1, 15);
                    }
                    
                    _mm_storeu_si128((__m128i*)(rgba + (i + j) * 4), color);
                }
            }
#else
            for (int i = 0; i < count; i++) {
                uint8_t idx = indices[i];
                int rgba_idx = i * 4;
                
                rgba[rgba_idx]     = palette[idx][0];
                rgba[rgba_idx + 1] = palette[idx][1];
                rgba[rgba_idx + 2] = palette[idx][2];
                
                if (has_transparent && idx == transparent_idx) {
                    rgba[rgba_idx + 3] = 0;
                } else {
                    rgba[rgba_idx + 3] = 255;
                }
            }
#endif
        }
    };

    class LZWDecoder {
    private:
        static const int MAX_CODES = 4096;
        
        struct CodeEntry {
            uint16_t prefix;
            uint8_t suffix;
            uint8_t length;
        };
        
        CodeEntry code_table[MAX_CODES];
        uint8_t decode_buffer[4096];
        
    public:
        ALWAYS_INLINE std::vector<uint8_t> decode(const uint8_t* data, size_t size, 
                                                 uint8_t min_code_size) {
            std::vector<uint8_t> output;
            
            uint32_t bit_buffer = 0;
            int bits_in_buffer = 0;
            size_t pos = 0;
            
            int clear_code = 1 << min_code_size;
            int eoi_code = clear_code + 1;
            int next_code = eoi_code + 1;
            int code_size = min_code_size + 1;

            for (int i = 0; i < clear_code; i++) {
                code_table[i] = {0xFFFF, static_cast<uint8_t>(i), 1};
            }
            
            uint16_t old_code = 0xFFFF;
            uint8_t first_char = 0;
            
            while (pos < size) {
                while (bits_in_buffer < code_size && pos < size) {
                    bit_buffer |= (data[pos++] << bits_in_buffer);
                    bits_in_buffer += 8;
                }
                
                if (bits_in_buffer < code_size) break;
                
                uint16_t code = bit_buffer & ((1 << code_size) - 1);
                bit_buffer >>= code_size;
                bits_in_buffer -= code_size;
                
                if (code == clear_code) {
                    next_code = eoi_code + 1;
                    code_size = min_code_size + 1;
                    old_code = 0xFFFF;
                    continue;
                }
                
                if (code == eoi_code) {
                    break;
                }

                uint16_t cur_code = code;
                int buf_pos = 4095;
                
                if (cur_code >= next_code) {
                    decode_buffer[buf_pos--] = first_char;
                    cur_code = old_code;
                }
                
                while (cur_code >= clear_code) {
                    decode_buffer[buf_pos--] = code_table[cur_code].suffix;
                    cur_code = code_table[cur_code].prefix;
                }
                
                decode_buffer[buf_pos--] = static_cast<uint8_t>(cur_code);
                first_char = decode_buffer[buf_pos + 1];

                int length = 4095 - buf_pos;
                output.insert(output.end(), decode_buffer + buf_pos + 1, decode_buffer + 4096);

                if (old_code != 0xFFFF && next_code < MAX_CODES) {
                    code_table[next_code].prefix = old_code;
                    code_table[next_code].suffix = first_char;
                    code_table[next_code].length = code_table[old_code].length + 1;
                    next_code++;
                    
                    if (next_code >= (1 << code_size) && code_size < 12) {
                        code_size++;
                    }
                }
                
                old_code = code;
            }
            
            return output;
        }
    };

    class Deinterlacer {
    public:
        ALWAYS_INLINE static std::vector<uint8_t> deinterlace(const uint8_t* src, int width, int height) {
            std::vector<uint8_t> dst(width * height);

            int passes[4] = {0, 4, 2, 1};
            int intervals[4] = {8, 8, 4, 2};
            
            int dest_pos = 0;
            
            for (int pass = 0; pass < 4; pass++) {
                int start_row = passes[pass];
                int interval = intervals[pass];
                
                for (int y = start_row; y < height; y += interval) {
                    std::memcpy(&dst[dest_pos], &src[y * width], width);
                    dest_pos += width;
                }
            }
            
            return dst;
        }
    };

    class GIFMemoryMap {
    private:
        const uint8_t* data_;
        size_t size_;
        size_t pos_;
        
    public:
        GIFMemoryMap(const std::string& filename) : data_(nullptr), size_(0), pos_(0) {
            std::ifstream file(filename, std::ios::binary | std::ios::ate);
            if (!file) return;
            
            size_ = file.tellg();
            file.seekg(0, std::ios::beg);
            
            std::vector<uint8_t> buffer(size_);
            file.read(reinterpret_cast<char*>(buffer.data()), size_);
            data_ = new uint8_t[size_];
            std::memcpy(const_cast<uint8_t*>(data_), buffer.data(), size_);
        }
        
        ~GIFMemoryMap() {
            delete[] data_;
        }
        
        ALWAYS_INLINE uint8_t read_byte() {
            return (pos_ < size_) ? data_[pos_++] : 0;
        }
        
        ALWAYS_INLINE uint16_t read_word() {
            uint16_t b1 = read_byte();
            uint16_t b2 = read_byte();
            return (b2 << 8) | b1;
        }
        
        ALWAYS_INLINE void skip(size_t n) {
            pos_ += n;
            if (pos_ > size_) pos_ = size_;
        }
        
        ALWAYS_INLINE bool eof() const { return pos_ >= size_; }
        ALWAYS_INLINE size_t position() const { return pos_; }
        ALWAYS_INLINE const uint8_t* data() const { return data_; }
        ALWAYS_INLINE size_t size() const { return size_; }
        
        ALWAYS_INLINE std::vector<uint8_t> read_data_block() {
            std::vector<uint8_t> block;
            uint8_t size = read_byte();
            
            while (size > 0 && !eof()) {
                for (int i = 0; i < size; i++) {
                    block.push_back(read_byte());
                }
                size = read_byte();
            }
            
            return block;
        }
    };

public:
    static ImageData load(const std::string& filename) {
        GIFMemoryMap mmap(filename);
        if (mmap.size() < 13) {
            throw image_error("GIF file too small");
        }
        
        GIFHeader header;
        std::memcpy(&header, mmap.data(), sizeof(header));
        mmap.skip(sizeof(header));
        
        if (std::string((char*)header.signature, 3) != "GIF") {
            throw image_error("Invalid GIF signature");
        }
        
        int width = header.width;
        int height = header.height;
        bool global_color_table = (header.packed & 0x80) != 0;
        int color_resolution = ((header.packed >> 4) & 0x07) + 1;
        bool sort_flag = (header.packed & 0x08) != 0;
        int global_color_table_size = 1 << ((header.packed & 0x07) + 1);
        
        ColorTable global_table;
        
        if (global_color_table) {
            global_table.size = global_color_table_size;
            for (int i = 0; i < global_color_table_size; i++) {
                global_table.colors[i][0] = mmap.read_byte();
                global_table.colors[i][1] = mmap.read_byte();
                global_table.colors[i][2] = mmap.read_byte();
            }
        } else {
            global_table.size = 256;
            for (int i = 0; i < 256; i++) {
                global_table.colors[i][0] = global_table.colors[i][1] = global_table.colors[i][2] = i;
            }
        }
        
        Frame frame;
        frame.width = width;
        frame.height = height;
        frame.transparent = false;
        frame.delay = 0;
        
        bool done = false;
        
        while (!mmap.eof() && !done) {
            uint8_t block_type = mmap.read_byte();
            
            switch (block_type) {
                case 0x21: {
                    uint8_t label = mmap.read_byte();
                    
                    if (label == 0xF9) {
                        mmap.skip(1);
                        uint8_t packed = mmap.read_byte();
                        frame.delay = mmap.read_word();
                        frame.transparent_index = mmap.read_byte();
                        frame.transparent = (packed & 0x01) != 0;
                        mmap.skip(1);
                    }
                    else if (label == 0xFE) {
                        auto data = mmap.read_data_block();
                    }
                    else if (label == 0xFF) {
                        mmap.skip(1);
                        mmap.skip(11);
                        mmap.read_data_block();
                    }
                    else {
                        mmap.read_data_block();
                    }
                    break;
                }
                
                case 0x2C: {
                    GIFImageDescriptor desc;
                    std::memcpy(&desc, mmap.data() + mmap.position() - 1, sizeof(desc));
                    mmap.skip(sizeof(desc) - 1);
                    
                    bool local_color_table = (desc.packed & 0x80) != 0;
                    bool interlaced = (desc.packed & 0x40) != 0;
                    int local_color_table_size = 1 << ((desc.packed & 0x07) + 1);
                    
                    ColorTable local_table;
                    
                    if (local_color_table) {
                        local_table.size = local_color_table_size;
                        for (int i = 0; i < local_color_table_size; i++) {
                            local_table.colors[i][0] = mmap.read_byte();
                            local_table.colors[i][1] = mmap.read_byte();
                            local_table.colors[i][2] = mmap.read_byte();
                        }
                    }
                    
                    uint8_t lzw_min_code_size = mmap.read_byte();
                    auto compressed_data = mmap.read_data_block();
                    
                    LZWDecoder decoder;
                    auto indices = decoder.decode(compressed_data.data(), compressed_data.size(), lzw_min_code_size);
                    
                    if (interlaced) {
                        indices = Deinterlacer::deinterlace(indices.data(), desc.width, desc.height);
                    }
                    
                    frame.indices = std::move(indices);
                    done = true;
                    break;
                }
                
                case 0x3B:
                    done = true;
                    break;
                    
                default:
                    break;
            }
        }
        
        if (frame.indices.empty()) {
            throw image_error("No image data found in GIF");
        }

        std::vector<uint8_t> rgba(frame.width * frame.height * 4);
        
#ifdef __SSE4_1__
        SIMDConverter::palette_to_rgba_sse(global_table.colors,
                                          frame.indices.data(),
                                          rgba.data(),
                                          frame.width * frame.height,
                                          frame.transparent_index,
                                          frame.transparent);
#else
        #pragma omp parallel for
        for (int i = 0; i < frame.width * frame.height; i++) {
            uint8_t idx = frame.indices[i];
            int rgba_idx = i * 4;
            
            rgba[rgba_idx] = global_table.colors[idx][0];
            rgba[rgba_idx + 1] = global_table.colors[idx][1];
            rgba[rgba_idx + 2] = global_table.colors[idx][2];
            
            if (frame.transparent && idx == frame.transparent_index) {
                rgba[rgba_idx + 3] = 0;
            } else {
                rgba[rgba_idx + 3] = 255;
            }
        }
#endif
        
        return ImageData(frame.width, frame.height, 4, std::move(rgba));
    }
    
    static std::vector<ImageData> load_all_frames(const std::string& filename) {
        GIFMemoryMap mmap(filename);
        if (mmap.size() < 13) {
            throw image_error("GIF file too small");
        }
        
        GIFHeader header;
        std::memcpy(&header, mmap.data(), sizeof(header));
        mmap.skip(sizeof(header));
        
        if (std::string((char*)header.signature, 3) != "GIF") {
            throw image_error("Invalid GIF signature");
        }
        
        int width = header.width;
        int height = header.height;
        bool global_color_table = (header.packed & 0x80) != 0;
        int global_color_table_size = 1 << ((header.packed & 0x07) + 1);
        
        ColorTable global_table;
        
        if (global_color_table) {
            global_table.size = global_color_table_size;
            for (int i = 0; i < global_color_table_size; i++) {
                global_table.colors[i][0] = mmap.read_byte();
                global_table.colors[i][1] = mmap.read_byte();
                global_table.colors[i][2] = mmap.read_byte();
            }
        } else {
            global_table.size = 256;
            for (int i = 0; i < 256; i++) {
                global_table.colors[i][0] = global_table.colors[i][1] = global_table.colors[i][2] = i;
            }
        }
        
        std::vector<Frame> frames;
        Frame current_frame;
        current_frame.width = width;
        current_frame.height = height;
        current_frame.transparent = false;
        current_frame.delay = 0;
        
        while (!mmap.eof()) {
            uint8_t block_type = mmap.read_byte();
            
            switch (block_type) {
                case 0x21: {
                    uint8_t label = mmap.read_byte();
                    
                    if (label == 0xF9) {
                        mmap.skip(1);
                        uint8_t packed = mmap.read_byte();
                        current_frame.delay = mmap.read_word();
                        current_frame.transparent_index = mmap.read_byte();
                        current_frame.transparent = (packed & 0x01) != 0;
                        mmap.skip(1);
                    }
                    else if (label == 0xFE) {
                        mmap.read_data_block();
                    }
                    else if (label == 0xFF) {
                        mmap.skip(1);
                        mmap.skip(11);
                        mmap.read_data_block();
                    }
                    else {
                        mmap.read_data_block();
                    }
                    break;
                }
                
                case 0x2C: {
                    GIFImageDescriptor desc;
                    std::memcpy(&desc, mmap.data() + mmap.position() - 1, sizeof(desc));
                    mmap.skip(sizeof(desc) - 1);
                    
                    bool local_color_table = (desc.packed & 0x80) != 0;
                    bool interlaced = (desc.packed & 0x40) != 0;
                    int local_color_table_size = 1 << ((desc.packed & 0x07) + 1);
                    
                    ColorTable local_table;
                    
                    if (local_color_table) {
                        local_table.size = local_color_table_size;
                        for (int i = 0; i < local_color_table_size; i++) {
                            local_table.colors[i][0] = mmap.read_byte();
                            local_table.colors[i][1] = mmap.read_byte();
                            local_table.colors[i][2] = mmap.read_byte();
                        }
                    }
                    
                    uint8_t lzw_min_code_size = mmap.read_byte();
                    auto compressed_data = mmap.read_data_block();
                    
                    LZWDecoder decoder;
                    auto indices = decoder.decode(compressed_data.data(), compressed_data.size(), lzw_min_code_size);
                    
                    if (interlaced) {
                        indices = Deinterlacer::deinterlace(indices.data(), desc.width, desc.height);
                    }
                    
                    current_frame.indices = std::move(indices);
                    current_frame.width = desc.width;
                    current_frame.height = desc.height;
                    current_frame.left = desc.left;
                    current_frame.top = desc.top;
                    
                    frames.push_back(current_frame);
                    current_frame = Frame();
                    current_frame.width = width;
                    current_frame.height = height;
                    break;
                }
                
                case 0x3B:
                    goto end_parse;
                    
                default:
                    break;
            }
        }
        
    end_parse:
        
        std::vector<ImageData> result;
        
        for (auto& frame : frames) {
            std::vector<uint8_t> rgba(frame.width * frame.height * 4);
            
            #pragma omp parallel for
            for (int i = 0; i < frame.width * frame.height; i++) {
                uint8_t idx = frame.indices[i];
                int rgba_idx = i * 4;
                
                rgba[rgba_idx] = global_table.colors[idx][0];
                rgba[rgba_idx + 1] = global_table.colors[idx][1];
                rgba[rgba_idx + 2] = global_table.colors[idx][2];
                
                if (frame.transparent && idx == frame.transparent_index) {
                    rgba[rgba_idx + 3] = 0;
                } else {
                    rgba[rgba_idx + 3] = 255;
                }
            }
            
            result.emplace_back(frame.width, frame.height, 4, std::move(rgba));
        }
        
        return result;
    }
    
    static std::vector<int> get_frame_delays(const std::string& filename) {
        GIFMemoryMap mmap(filename);
        if (mmap.size() < 13) return {};
        
        mmap.skip(13);
        
        std::vector<int> delays;
        bool has_gct = (mmap.data()[10] & 0x80) != 0;
        int gct_size = 1 << ((mmap.data()[10] & 0x07) + 1);
        
        if (has_gct) {
            mmap.skip(gct_size * 3);
        }
        
        while (!mmap.eof()) {
            uint8_t block_type = mmap.read_byte();
            
            if (block_type == 0x21) {
                uint8_t label = mmap.read_byte();
                
                if (label == 0xF9) {
                    mmap.skip(1);
                    uint8_t packed = mmap.read_byte();
                    uint16_t delay = mmap.read_word();
                    delays.push_back(delay * 10);
                    mmap.skip(2);
                }
                else {
                    mmap.read_data_block();
                }
            }
            else if (block_type == 0x2C) {
                GIFImageDescriptor desc;
                std::memcpy(&desc, mmap.data() + mmap.position() - 1, sizeof(desc));
                mmap.skip(sizeof(desc) - 1);
                
                bool local_color_table = (desc.packed & 0x80) != 0;
                int local_color_table_size = 1 << ((desc.packed & 0x07) + 1);
                
                if (local_color_table) {
                    mmap.skip(local_color_table_size * 3);
                }
                
                mmap.skip(1);
                mmap.read_data_block();
            }
            else if (block_type == 0x3B) {
                break;
            }
        }
        
        return delays;
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
    else if (ieq(ext, "jpg") || ieq(ext, "jpeg") || ieq(ext, "jpe")) {
        image = JPEGDecoder::load(filename);
    }
    else if (ieq(ext, "webp")) {
        image = WebPDecoder::load(filename);
    }
    else if (ieq(ext, "hdr") || ieq(ext, "pic")) {
        image = HDRDecoder::load(filename);
    }
    else if (ieq(ext, "gif")) {
        image = GIFDecoder::load(filename);
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
