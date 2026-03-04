#include <merkle.h>
#include <NTL/ZZ.h>
#include <cstring>
#include <cassert>

using namespace NTL;

// ============================================================
// FNV-1a hash stretched to 32 bytes (prototype only).
// For production: replace with OpenSSL SHA-256.
//   #include <openssl/sha.h>
//   SHA256(data, len, out);
// ============================================================
static void simple_hash_256(const unsigned char* data, size_t len, unsigned char out[32]) {
    uint64_t h1 = 14695981039346656037ULL;
    uint64_t h2 = 14695981039346656037ULL ^ 0xDEADBEEF;
    uint64_t h3 = 14695981039346656037ULL ^ 0xCAFEBABE;
    uint64_t h4 = 14695981039346656037ULL ^ 0x12345678;
    for (size_t i = 0; i < len; i++) {
        h1 ^= data[i]; h1 *= 1099511628211ULL;
        h2 ^= data[i]; h2 *= 1099511628213ULL;
        h3 ^= data[i]; h3 *= 1099511628249ULL;
        h4 ^= data[i]; h4 *= 1099511628307ULL;
    }
    memcpy(out,      &h1, 8);
    memcpy(out + 8,  &h2, 8);
    memcpy(out + 16, &h3, 8);
    memcpy(out + 24, &h4, 8);
}

std::vector<unsigned char> sha256(const unsigned char* data, size_t len) {
    std::vector<unsigned char> out(32);
    simple_hash_256(data, len, out.data());
    return out;
}

std::vector<unsigned char> sha256(const std::vector<unsigned char>& data) {
    return sha256(data.data(), data.size());
}

// ============================================================
// FIX: Renamed local variable from 'coeff' to 'c_zz'
// to avoid shadowing NTL::coeff() free function.
// ============================================================
std::vector<unsigned char> serialize_ZZ_pE(const ZZ_pE& elem) {
    long d = ZZ_pE::degree();
    const ZZ_pX& poly = rep(elem);
    std::vector<unsigned char> bytes;
    long num_bytes = NumBytes(ZZ_p::modulus());

    for (long i = 0; i < d; i++) {
        ZZ_p ci;
        if (!IsZero(poly) && deg(poly) >= i) {
            ci = coeff(poly, i);  // NTL::coeff() — no longer shadowed
        } else {
            clear(ci);
        }
        ZZ c_zz = rep(ci);       // 'c_zz' instead of 'coeff'
        std::vector<unsigned char> buf(num_bytes, 0);
        BytesFromZZ(buf.data(), c_zz, num_bytes);
        bytes.insert(bytes.end(), buf.begin(), buf.end());
    }
    return bytes;
}

std::vector<unsigned char> hash_column(const ZZ_pE* column, long len) {
    std::vector<unsigned char> data;
    for (long i = 0; i < len; i++) {
        auto s = serialize_ZZ_pE(column[i]);
        data.insert(data.end(), s.begin(), s.end());
    }
    return sha256(data);
}

std::vector<unsigned char> hash_pair(const std::vector<unsigned char>& left,
                                     const std::vector<unsigned char>& right) {
    std::vector<unsigned char> combined;
    combined.insert(combined.end(), left.begin(), left.end());
    combined.insert(combined.end(), right.begin(), right.end());
    return sha256(combined);
}

MerkleTree build_merkle_tree(const std::vector<std::vector<unsigned char>>& leaf_hashes) {
    MerkleTree tree;
    tree.num_leaves = leaf_hashes.size();

    if (tree.num_leaves <= 1) {
        tree.depth = 0;
        tree.hashes = leaf_hashes;
        tree.root = leaf_hashes.empty()
            ? std::vector<unsigned char>(32, 0)
            : leaf_hashes[0];
        return tree;
    }

    // Pad to next power of 2
    long p2 = 1;
    tree.depth = 0;
    while (p2 < tree.num_leaves) { p2 <<= 1; tree.depth++; }

    long total = (2 * p2) - 1;
    tree.hashes.resize(total, std::vector<unsigned char>(32, 0));

    // Fill leaves
    for (long i = 0; i < tree.num_leaves; i++) {
        tree.hashes[i] = leaf_hashes[i];
    }

    // Build tree bottom-up
    long offset = 0;
    for (long width = p2; width > 1; width >>= 1) {
        for (long i = 0; i < width; i += 2) {
            tree.hashes[offset + width + i / 2] =
                hash_pair(tree.hashes[offset + i], tree.hashes[offset + i + 1]);
        }
        offset += width;
    }

    tree.root = tree.hashes.back();
    return tree;
}

std::vector<std::vector<unsigned char>> get_merkle_path(const MerkleTree& tree, long leaf_idx) {
    std::vector<std::vector<unsigned char>> path;
    long p2 = 1;
    while (p2 < tree.num_leaves) p2 <<= 1;

    long offset = 0;
    long idx = leaf_idx;
    for (long width = p2; width > 1; width >>= 1) {
        long neighbor = idx ^ 1;
        path.push_back(tree.hashes[offset + neighbor]);
        offset += width;
        idx >>= 1;
    }
    return path;
}

bool verify_merkle_path(const std::vector<unsigned char>& leaf_hash,
                        long leaf_idx,
                        const std::vector<std::vector<unsigned char>>& path,
                        const std::vector<unsigned char>& root) {
    auto current = leaf_hash;
    long idx = leaf_idx;
    for (size_t level = 0; level < path.size(); level++) {
        if ((idx & 1) == 0) {
            current = hash_pair(current, path[level]);
        } else {
            current = hash_pair(path[level], current);
        }
        idx >>= 1;
    }
    return current == root;
}