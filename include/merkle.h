#ifndef MERKLE_H
#define MERKLE_H

#include <NTL/ZZ_pE.h>
#include <vector>
#include <string>

using namespace NTL;

// ============================================================
// SHA-256 based Merkle tree for ZZ_pE columns
// ============================================================

// Serialize a ZZ_pE element to bytes (deterministic)
std::vector<unsigned char> serialize_ZZ_pE(const ZZ_pE& elem);

// SHA-256 hash of raw bytes
std::vector<unsigned char> sha256(const unsigned char* data, size_t len);
std::vector<unsigned char> sha256(const std::vector<unsigned char>& data);

// Hash a vector of ZZ_pE elements (one column of the encoded matrix)
std::vector<unsigned char> hash_column(const ZZ_pE* column, long len);

// Hash two children to produce parent
std::vector<unsigned char> hash_pair(const std::vector<unsigned char>& left,
                                     const std::vector<unsigned char>& right);

// H_M: Hash M hashes together (M-ary hash function for distributed protocol)
// h[k] = H_M(h^(0)[k], h^(1)[k], ..., h^(M-1)[k])
std::vector<unsigned char> hash_M(const std::vector<std::vector<unsigned char>>& hashes);

// Build a Merkle tree from leaf hashes
// Returns: flat array of hashes. Layout: [leaves..., level1..., level2..., ..., root]
// Total size: (2 << depth) - 1, where depth = ceil(log2(num_leaves))
struct MerkleTree {
    std::vector<std::vector<unsigned char>> hashes;
    long num_leaves;
    long depth;
    std::vector<unsigned char> root;
};

MerkleTree build_merkle_tree(const std::vector<std::vector<unsigned char>>& leaf_hashes);

// Get Merkle path for a given leaf index
std::vector<std::vector<unsigned char>> get_merkle_path(const MerkleTree& tree, long leaf_idx);

// Verify a Merkle path
bool verify_merkle_path(const std::vector<unsigned char>& leaf_hash,
                        long leaf_idx,
                        const std::vector<std::vector<unsigned char>>& path,
                        const std::vector<unsigned char>& root);

#endif