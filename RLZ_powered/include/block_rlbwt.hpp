#pragma once

#include <cmath>
#include <random>
#include <cassert>
#include <cstring>
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <cstdint>
#include <bits/fs_fwd.h>
#include <bitset>

namespace bbwt {
template <class bwt_type>
class block_rlbwt_builder {
   private:
    static const constexpr uint64_t BLOCKS_IN_SUPER_BLOCK =
        bwt_type::super_block_type::blocks;

    typedef typename bwt_type::alphabet_type alphabet_type;
    typedef typename bwt_type::block_alphabet_type block_alphabet_type;
    typedef typename bwt_type::block_type block_type;

    uint64_t char_counts_[257];
    uint32_t run_count_;
    uint32_t dense_blocks_;
    std::string prefix_;
    std::string suffix_;
    std::vector<alphabet_type> block_counts_;
    std::vector<bool> block_reprs_;
    alphabet_type super_block_cumulative_;
    std::vector<uint64_t> block_offsets_;
    std::vector<uint8_t> sa_;
    block_alphabet_type block_cumulative_;
    uint64_t super_block_bytes_;
    uint64_t super_block_size_;
    uint64_t elems_;
    uint8_t* current_super_block_;
    uint8_t** scratch_;
    block_type current_block_;
    uint32_t block_elems_;
    uint32_t block_bytes_;
    uint32_t blocks_in_super_block_;
    std::fstream out_;

   public:
    block_rlbwt_builder(std::string out_file, std::string sa_file_name = "")
        : char_counts_(),
          run_count_(0),
          dense_blocks_(0),
          block_counts_(),
          block_reprs_(),
          super_block_cumulative_(),
          block_offsets_(),
          sa_(),
          block_cumulative_(),
          super_block_bytes_(sizeof(block_alphabet_type)),
          super_block_size_(
              BLOCKS_IN_SUPER_BLOCK *
              (block_type::min_size + sizeof(block_alphabet_type))),
          elems_(0),
          current_block_(),
          block_elems_(0),
          block_bytes_(0),
          blocks_in_super_block_(0) {
        size_t loc = out_file.find_last_of('.');
        if (loc == std::string::npos) {
            prefix_ = out_file;
            suffix_ = "";
        } else {
            prefix_ = out_file.substr(0, loc);
            suffix_ = out_file.substr(loc);
        }
        out_.open(prefix_ + "_data" + suffix_, std::ios::binary | std::ios::out);
        block_counts_.push_back(super_block_cumulative_);
        current_super_block_ = (uint8_t*)calloc(super_block_size_, 1);
        scratch_ =
            (uint8_t**)malloc(block_type::scratch_blocks * sizeof(uint8_t*));
        for (size_t i = 0; i < block_type::scratch_blocks; i++) {
            scratch_[i] = (uint8_t*)calloc(block_type::scratch_size(i), 1);
        }

        std::ifstream sa_file(sa_file_name, std::ios::binary);
        if (!sa_file) {
            throw std::runtime_error("Cannot open sa file: " + sa_file_name);
        }

        sa_file.seekg(0, std::ios::end);
        std::size_t sa_file_size = sa_file.tellg();
        sa_file.seekg(0, std::ios::beg);

        if (sa_file_size % sizeof(uint8_t) != 0) {
            throw std::runtime_error("File size is not a multiple of uint64_t.");
        }

        sa_.resize(sa_file_size / sizeof(uint8_t));
        sa_file.read(reinterpret_cast<char*>(sa_.data()), sa_file_size);

        if (!sa_file) {
            throw std::runtime_error("Error while reading file.");
        }
    }

    void append(uint8_t head, uint32_t length) {
        char_counts_[head] += length;
        head = alphabet_type::convert(head);
        while (length) {
            run_count_++;
            if (length + block_elems_ < bwt_type::cap) {
                block_cumulative_.add(head, length);
                super_block_cumulative_.add(head, length);
                block_bytes_ = current_block_.append(head, length, scratch_);
                block_elems_ += length;
                elems_ += length;
                return;
            } else if (length + block_elems_ == bwt_type::cap) [[unlikely]] {
                block_cumulative_.add(head, length);
                super_block_cumulative_.add(head, length);
                block_bytes_ = current_block_.append(head, length, scratch_);
                elems_ += length;
                commit();
                return;
            } else {
                uint32_t fill = bwt_type::cap - block_elems_;
                block_cumulative_.add(head, fill);
                super_block_cumulative_.add(head, fill);
                block_bytes_ = current_block_.append(head, fill, scratch_);
                elems_ += fill;
                commit();
                length -= fill;
            }
        }
    }

    void finalize() {
        if (block_elems_) {
            commit(true);
        }
        if (blocks_in_super_block_) {
            write_super_block();
        }
        //out_.close();
        uint64_t p_v = 0;
        for (size_t i = 0; i < 257; i++) {
            uint64_t tmp = char_counts_[i];
            char_counts_[i] = p_v;
            p_v += tmp;
        }
        write_root();
    }

    void gen_queries(std::ostream& out, uint32_t n_queries) {
        uint64_t char_count = char_counts_[256];
        std::vector<uint8_t> chars;
        for (uint16_t i = 0; i < 256; i++) {
            if (char_counts_[i + 1] > char_counts_[i]) {
                chars.push_back(uint8_t(i));
            }
        }

        std::mt19937 mt;
        std::uniform_int_distribution<unsigned long long> i_gen(0, char_count - 1);
        std::uniform_int_distribution<uint8_t> c_gen(0, chars.size() - 1);

        for (uint32_t i = 0; i < n_queries; i++) {
            uint64_t idx = i_gen(mt);
            uint8_t c = chars[c_gen(mt)];
            bool dense = block_reprs_[idx / bwt_type::cap];
            out.write(reinterpret_cast<char*>(&idx), sizeof(uint64_t));
            out.write(reinterpret_cast<char*>(&c), sizeof(uint8_t));
            out.write(reinterpret_cast<char*>(&dense), sizeof(bool));
        }
    }

   private:
    void write_super_block() {
        uint64_t file_bytes =
            super_block_bytes_ +
            sizeof(uint64_t) * bwt_type::super_block_type::blocks;
        out_.write(reinterpret_cast<char*>(&file_bytes), sizeof(uint64_t));

        uint64_t block_offsets =
            sizeof(uint64_t) * block_offsets_.size();
        out_.write(reinterpret_cast<char*>(&block_offsets), sizeof(uint64_t));
        out_.write(reinterpret_cast<char*>(block_offsets_.data()),
                  sizeof(uint64_t) * block_offsets_.size());
        uint64_t padding_bytes = sizeof(uint64_t) * ( bwt_type::super_block_type::blocks - block_offsets_.size());
        out_.write(reinterpret_cast<char*>(&padding_bytes), sizeof(uint64_t));
        out_.write(reinterpret_cast<char*>(&super_block_bytes_), sizeof(uint64_t));
        out_.write(reinterpret_cast<char*>(current_super_block_),
                  super_block_bytes_);
        block_counts_.push_back(super_block_cumulative_);

        block_cumulative_.clear();
        block_offsets_.clear();
        std::memset(current_super_block_, 0,
                    sizeof(uint8_t) * super_block_size_);
        super_block_bytes_ = sizeof(block_alphabet_type);
        blocks_in_super_block_ = 0;
    }

    void commit(bool last_block = false) {
        if (run_count_ >= 4 * std::log2(bwt_type::cap)) {
            block_reprs_.push_back(true);
            dense_blocks_++;
        } else {
            block_reprs_.push_back(false);
        }
        run_count_ = 0;
        if (super_block_bytes_ + block_bytes_ + block_alphabet_type::size() > super_block_size_) {
            uint64_t new_size = super_block_bytes_ + block_bytes_;
            if (!last_block) {
                new_size +=
                    (BLOCKS_IN_SUPER_BLOCK - blocks_in_super_block_ - 1) *
                    (block_type::min_size + sizeof(block_alphabet_type));
            }
            current_super_block_ =
                (uint8_t*)realloc(current_super_block_, new_size);
            super_block_size_ = new_size;
        }
        
        block_offsets_.push_back(super_block_bytes_);
        if constexpr (block_type::has_members) {
            std::memcpy(current_super_block_ + super_block_bytes_, &current_block_,
                    sizeof(block_type));
        }
        block_type* b = reinterpret_cast<block_type*>(current_super_block_ +
                                                      super_block_bytes_);
        super_block_bytes_ += b->commit(scratch_);

        for (size_t i = 0; i < block_type::scratch_blocks; i++) {
            std::memset(scratch_[i], 0, block_type::scratch_size(i));
        }
        current_block_.clear();
        block_elems_ = 0;
        block_bytes_ = 0;
        blocks_in_super_block_++;

        if (!last_block && blocks_in_super_block_ < BLOCKS_IN_SUPER_BLOCK)
            [[likely]] {
            std::memcpy(current_super_block_ + super_block_bytes_,
                        &block_cumulative_, block_alphabet_type::size());
            super_block_bytes_ += block_alphabet_type::size();
        } else {
            write_super_block();
        }
    }

    void write_root() {
        uint64_t n_blocks = block_counts_.size() - 1;

        std::cerr << "Writing \"root\" of " << block_type::cap << "-sb-rlbwt to file\n"
                  << " Seen " << n_blocks << " super blocks\n"
                  << " containing a total of " << elems_ << " elements\n"
                  << " " << dense_blocks_ << " dense blocks out of "
                  << block_reprs_.size() << " total blocks ("
                  << 100.0 * dense_blocks_ / block_reprs_.size() << " %)"
                  << std::endl;

        std::fstream out;
        out.open(prefix_ + suffix_, std::ios::binary | std::ios::out);
        alphabet_type::write_statics(out);
        block_alphabet_type::write_statics(out);
        bwt_type::super_block_type::write_statics(out);
        block_type::write_statics(out);
        uint64_t bytes = sizeof(alphabet_type) * block_counts_.size();
        out.write(reinterpret_cast<char*>(&bytes), sizeof(uint64_t));
        out.write(reinterpret_cast<char*>(&elems_), sizeof(uint64_t));
        out.write(reinterpret_cast<char*>(&n_blocks), sizeof(uint64_t));
        out.write(reinterpret_cast<char*>(block_counts_.data()), bytes);
        out.write(reinterpret_cast<char*>(char_counts_), sizeof(uint64_t) * 257);
        uint64_t sa_size = sizeof(uint8_t) * sa_.size();
        out.write(reinterpret_cast<char*>(&sa_size), sizeof(uint64_t));
        out.write(reinterpret_cast<char*>(sa_.data()), sa_size);
        out.close();
    }
};

template <class super_block_type_, class alphabet_type_>
class block_rlbwt {
   public:
    typedef super_block_type_ super_block_type;
    typedef alphabet_type_ alphabet_type;
    typedef block_rlbwt_builder<block_rlbwt> builder;

   private:
    static const constexpr uint64_t SUPER_BLOCK_ELEMS = uint64_t(1) << 32;
    uint64_t size_;
    uint64_t block_count_;
    uint64_t bytes_;
    uint64_t char_counts_[257];
    uint8_t* p_sums_;
    std::vector<super_block_type*> s_blocks_;

   public:
    static const constexpr uint32_t cap = super_block_type::cap;

    std::vector<uint64_t> sa_;

    static_assert(SUPER_BLOCK_ELEMS % cap == 0);

    typedef typename super_block_type_::block_type block_type;
    typedef typename block_type::alphabet_type block_alphabet_type;

    block_rlbwt(std::string path) : bytes_(sizeof(block_rlbwt)) {
        std::fstream in_file;
        in_file.open(path, std::ios::binary | std::ios::in);
        if (in_file.fail()) {
            std::cerr << " -> Failed" << std::endl;
            exit(1);
        }
        bytes_ += alphabet_type::load_statics(in_file);
        bytes_ += block_alphabet_type::load_statics(in_file);
        bytes_ += super_block_type::load_statics(in_file);
        bytes_ += block_type::load_statics(in_file);
        uint64_t data_bytes;
        in_file.read(reinterpret_cast<char*>(&data_bytes), sizeof(uint64_t));
        in_file.read(reinterpret_cast<char*>(&size_), sizeof(uint64_t));
        in_file.read(reinterpret_cast<char*>(&block_count_), sizeof(uint64_t));
#ifdef VERB
        std::cerr << data_bytes << " bytes of data\n"
                  << size_ << " logical elements\n"
                  << "in " << block_count_ << " super blocks" << std::endl;
#endif
        p_sums_ = (uint8_t*)std::malloc(data_bytes);
        in_file.read(reinterpret_cast<char*>(p_sums_), data_bytes);
        bytes_ += data_bytes;
        in_file.read(reinterpret_cast<char*>(char_counts_), sizeof(uint64_t) * 257);
        uint64_t sa_size;
        in_file.read(reinterpret_cast<char *>(&sa_size), sizeof(uint64_t));
        sa_size /= sizeof(uint8_t);
        std::vector<uint8_t> tmp_sa(sa_size);
        in_file.read(reinterpret_cast<char *>(tmp_sa.data()), sa_size);
        in_file.close();

        uint8_t data_size = (sa_size/size_);
        sa_.resize(size_);

        size_t index = 0;
        while (index < sa_size) {
            uint64_t conjugate = 0;
            for (uint8_t j = 0; j < data_size; j++) {
                conjugate += tmp_sa[index+j] << (8*j);
            }
            sa_[index/data_size] = conjugate;
            index += data_size;
        }

        std::string prefix;
        std::string suffix;
        size_t loc = path.find_last_of('.');
        if (loc == std::string::npos) {
            prefix = path;
            suffix = "";
        } else {
            prefix = path.substr(0, loc);
            suffix = path.substr(loc);
        }
        in_file.open(prefix + "_data" + suffix, std::ios::binary | std::ios::in);
        if (in_file.fail()) {
            std::cerr << " -> Failed" << std::endl;
            exit(1);
        }
        for (uint64_t i = 1; i <= block_count_; i++) {
            s_blocks_.push_back(read_super_block(in_file));
        }
        bytes_ += s_blocks_.size() * sizeof(super_block_type*);
        in_file.close();
    }

    block_rlbwt() = delete;
    block_rlbwt(const block_rlbwt& other) = delete;
    block_rlbwt& operator=(const block_rlbwt& other) = delete;

    block_rlbwt(block_rlbwt&& other) {
        size_ = std::exchange(other.size_, 0);
        bytes_ = std::exchange(other.bytes_, 0);
        block_count_ = std::exchange(other.block_count_, 0);
        p_sums_ = std::exchange(other.p_sums_, nullptr);
        s_blocks_ =
            std::exchange(other.s_blocks_, std::vector<super_block_type*>());
        std::memcpy(char_counts_, other.char_counts_, sizeof(uint64_t) * 257);
    }

    block_rlbwt& operator=(block_rlbwt&& other) {
        size_ = std::exchange(other.size_, 0);
        block_count_ = std::exchange(other.block_count_, 0);
        p_sums_ = std::exchange(other.p_sums_, nullptr);
        s_blocks_ =
            std::exchange(other.s_blocks_, std::vector<super_block_type*>());
        std::memcpy(char_counts_, other.char_counts_, sizeof(uint64_t) * 257);
        return *this;
    }

    ~block_rlbwt() {
        for (uint64_t i = 0; i < block_count_; i++) {
            std::free(s_blocks_[i]);
        }
        std::free(p_sums_);
    }

    uint8_t at(uint64_t i) const {
        if (i >= size_) [[unlikely]] {
            return 0;
        }
        uint64_t s_block_i = i / SUPER_BLOCK_ELEMS;

        return alphabet_type::revert(s_blocks_[s_block_i]->at(i % SUPER_BLOCK_ELEMS));
    }

    uint64_t count(const std::string_view& pattern, uint8_t offset, uint8_t code_size) const {
        uint8_t c = 0;

        for (uint8_t i = offset; i > 1; i--) {
            char ch = pattern[pattern.size() -i];
            c += ((ch >= 'C') + (ch >= 'G') + (ch >= 'T'));
            c <<= 2;
        }
        if (code_size > 1) {
            char ch = pattern[pattern.size() -1];
            c += ((ch >= 'C') + (ch >= 'G') + (ch >= 'T'));
            c <<= (code_size-offset)*2;
        } else {
            c = pattern[pattern.size() -1];
        }

        uint64_t a = char_counts_[c];
        uint64_t b = char_counts_[c + (1ULL<<((code_size-offset)*2))];

        for (size_t i = pattern.size() - offset -1; i < pattern.size() && b > a; i -= code_size) {

            if (code_size > 1) {
                char ch = pattern[i];
                c = ((ch >= 'C') + (ch >= 'G') + (ch >= 'T'));
                for (uint8_t j = 1; j < code_size; j++) {
                    ch = pattern[i-j];
                    c += ((ch >= 'C') + (ch >= 'G') + (ch >= 'T'))<<(j*2);
                }
            } else {
                c = pattern[i];
            }

            a = rank(a, c);
            b = rank(b, c);

            if (b <= a) [[unlikely]] {
                return 0;
            }

            a += char_counts_[c];
            b += char_counts_[c];
        }

        return b-a;
    }

    /** Searching for the largest extension of a metacharacter in the bwt interval.
     * When searching for the largest matched prefix, the extension step may only cover a part of the metacharacter.
     * This function provides means for such extension.
     *
     * @param a start of the interval in question
     * @param b end of the interval in question (left-open interval)
     * @param c metacharacter
     * @return tuple (l, g), where l is the extended length of the matched suffix and g is its starting position
     */
    std::pair<uint64_t, uint64_t> largest_extension (const uint64_t a, const uint64_t b, uint8_t c) const {
        std::pair<uint64_t, uint64_t> res;

        uint64_t s_block_a = a / SUPER_BLOCK_ELEMS;
        uint64_t s_block_b = b / SUPER_BLOCK_ELEMS;

        if (s_block_a == s_block_b) {
            res = s_blocks_[s_block_a]->extension(a % SUPER_BLOCK_ELEMS, b % SUPER_BLOCK_ELEMS, c);
            return {res.first, (sa_[(res.second + s_block_a*SUPER_BLOCK_ELEMS)] + size_ - res.first) % size_};
        }


        res = s_blocks_[s_block_a]->extension(a % SUPER_BLOCK_ELEMS, SUPER_BLOCK_ELEMS-1, c);
        res = {res.first, res.second + s_block_a*SUPER_BLOCK_ELEMS};
        s_block_a++;
        while (s_block_a < s_block_b) {
            auto r = s_blocks_[s_block_a]->extension(0, SUPER_BLOCK_ELEMS-1, c);

            if (r.first > res.first) res = {r.first, r.second+s_block_a*SUPER_BLOCK_ELEMS};
            s_block_a++;
        }

        auto r = s_blocks_[s_block_b]->extension(0, b % SUPER_BLOCK_ELEMS, c);
        if (r.first > res.first) res = {r.first, r.second+s_block_b*SUPER_BLOCK_ELEMS};

        return {res.first, (sa_[res.second] + size_ - res.first) % size_};
    }

    /** Function to compute the parses - tuples (length, position) - of the input sequence against this reference
     * The parses are computed by greedy search of the largest matching suffix in a loop
     *
     * @param seq sequence to be parsed against this FM-index
     * @param parses vector of the resulting pairs (length, starting position)
     * @param code_size number of characters in the metacharacter
     */
    void parse_tuples (const std::string_view& seq, std::vector<std::tuple<uint64_t, uint64_t>>& parses, uint8_t code_size) {
        auto constexpr encode_base = [](unsigned char ch) noexcept ->
                            uint8_t { return static_cast<uint8_t>
                                ( (ch > 'A') + (ch > 'C') + (ch > 'G')); };
        uint8_t c = 0;
        std::pair<uint64_t, uint64_t> tuple;

        // case if sequence is too short
        if (seq.size() <= code_size) {
            for (uint8_t j = 0; j < seq.size(); j++) {
                c <<= 2;
                c += encode_base(seq[j]);
            }
            c <<= ((code_size-seq.size())*2);

            while (char_counts_[c+1] == char_counts_[c]) c++;
            parses.emplace_back(seq.size(), sa_[char_counts_[c]]);
            return;
        }

        uint64_t length = code_size;
        size_t i = seq.size() - code_size;

        if (code_size == 1) {
            c = seq[i];
        } else {
            c = (encode_base(seq[i])<<6) | (encode_base(seq[i+1])<<4) | (encode_base(seq[i+2])<<2) | (encode_base(seq[i+3]));
        }

        uint64_t a = char_counts_[c];
        uint64_t b = char_counts_[c + 1];

        uint64_t a_tmp = 0;
        uint64_t b_tmp = 0;

        for (; i >= code_size; i -= code_size) {
            //load the next character / metacharacter
            if (code_size == 1)
                c = seq[i-1];
            else {
                c = encode_base(seq[i-4])<<6 | encode_base(seq[i-3])<<4 | encode_base(seq[i-2])<<2 | encode_base(seq[i-1]);
            }

            // compute the rank of this metacharacter
            a_tmp = rank(a, c);
            b_tmp = rank(b, c);

            // check if the suffix can be extended by a whole metacharacter, else extend by a partial metacharacter
            if (b_tmp <= a_tmp) {
                // check if the suffix can be extended by a partial metacharacter
                if (code_size > 1) {
                    tuple = largest_extension(a, b-1, c);
                    parses.emplace_back(length+tuple.first, tuple.second);

                    // if it was extended, recompute the current position and metacharacter in the sequence
                    if (tuple.first > 0) {
                        i -= tuple.first;

                        if (i < code_size || i >= seq.size()) {
                            length = 0;
                            break;
                        }

                        c = encode_base(seq[i-4])<<6 | encode_base(seq[i-3])<<4 | encode_base(seq[i-2])<<2 | encode_base(seq[i-1]);
                    }

                } else {
                    parses.emplace_back(length, sa_[a]);
                }

                // set the correct start interval
                length = code_size;
                a = char_counts_[c];
                b = char_counts_[c+1];
                continue;
            }

            // continue with extension
            a = a_tmp + char_counts_[c];
            b = a - a_tmp + b_tmp;

            //update length
            length += code_size;
        }

        // check if the loop covered the whole sequence or there is a partial metacharacter left to cover
        if (i > 0) {
            c = 0;
            for (uint8_t j = 0; j < i; j++) {
                c += encode_base(seq[i -j -1])<<(j<<1);
            }

            // extend the suffix
            if (length > 0) {
                tuple = largest_extension(a, b-1, c);
                if (tuple.first >= i) {
                    parses.emplace_back(length+i, (tuple.second+(tuple.first-i)) % size_);
                    return;
                }
                parses.emplace_back(length, sa_[a]);
            }
            c <<= ((code_size-i)<<1);

            while (char_counts_[c+1] == char_counts_[c]) c++;
            parses.emplace_back(i, sa_[char_counts_[c]]);
        } else {
            parses.emplace_back(length, sa_[a]);
        }
    }


    /** Function to compute the parses - triples (length, position, character) - of the input sequence against this reference
     * The parses are computed by greedy search of the largest matching suffix in a loop
     *
     * @param seq sequence to be parsed against this FM-index
     * @param parses vector of the resulting pairs (length, starting position)
     * @param code_size number of characters in the metacharacter
     */
    void parse_triples (const std::string_view& seq, std::vector<std::tuple<uint64_t, uint64_t, uint8_t>>& parses, uint8_t code_size) {
        auto constexpr encode_base = [](unsigned char ch) noexcept ->
                            uint8_t { return static_cast<uint8_t>
                                ( (ch > 'A') + (ch > 'C') + (ch > 'G')); };

        uint8_t c = 0;
        std::pair<uint64_t, uint64_t> tuple;
        uint8_t last = seq[seq.size()-1];

        if (seq.size() <= code_size) {
            for (uint8_t j = 0; j < seq.size()-1; j++) {
                c <<= 2;
                c += encode_base(seq[j]);
            }
            c <<= ((code_size-seq.size()+1)<<1);

            while (char_counts_[c+1] == char_counts_[c]) c++;
            parses.emplace_back(seq.size()-1, sa_[char_counts_[c]], last);
            return;
        }


        uint64_t length = code_size;
        size_t i = seq.size() - code_size -1;

        if (code_size == 1) {
            c = seq[i];
        } else {
            c = encode_base(seq[i])<<6 | encode_base(seq[i+1])<<4 | encode_base(seq[i+2])<<2 | encode_base(seq[i+3]);
        }

        uint64_t a = char_counts_[c];
        uint64_t b = char_counts_[c + 1];

        uint64_t a_tmp = 0;
        uint64_t b_tmp = 0;

        for (; i >= code_size; i -= code_size) {

            //load the next metacharacter
            if (code_size == 1)
                c = seq[i-1];
            else {
                c = encode_base(seq[i-4])<<6 | encode_base(seq[i-3])<<4 | encode_base(seq[i-2])<<2 | encode_base(seq[i-1]);
            }

            // compute the rank of this metacharacter
            a_tmp = rank(a, c);
            b_tmp = rank(b, c);

            // check if the prefix can be extended by a whole metacharacter
            if (b_tmp <= a_tmp) {
                // check if the prefix can be extended by a partial metacharacter
                if (code_size > 1) {
                    tuple = largest_extension(a, b, c);
                    parses.emplace_back(tuple.first+length, tuple.second, last);

                    // recompute the current position and metacharacter in the sequence
                    i -= tuple.first+1;
                    last = seq[i];

                    if (i < code_size || i >= seq.size()) {
                        length = 0;
                        break;
                    }

                    c = encode_base(seq[i-4])<<6 | encode_base(seq[i-3])<<4 | encode_base(seq[i-2])<<2 | encode_base(seq[i-1]);
                } else {
                    parses.emplace_back(length, sa_[a], last);
                    i--;
                    c = seq[i-1];
                    last = seq[i];

                    if (i < code_size || i >= seq.size()) {
                        length = 0;
                        break;
                    }
                }


                // set the correct start interval
                length = code_size;
                a = char_counts_[c];
                b = char_counts_[c+1];
                continue;
            }

            // continue with extension
            a = a_tmp + char_counts_[c];
            b = a - a_tmp + b_tmp;

            //update length
            length += code_size;
        }

        // check if the loop covered the whole sequence or there is a partial metacharacter left to cover
        if (i > 0) {
            c = 0;
            for (uint8_t j = 0; j < i; j++) {
                c += encode_base(seq[i -j -1])<<(j<<1);
            }
            if (length > 0) {
                tuple = largest_extension(a, b-1, c);
                if (tuple.first >= i) {
                    parses.emplace_back(length+i, (tuple.second+(tuple.first-i)) % size_, last);
                    return;
                }
                parses.emplace_back(length, sa_[a], last);
            }
            c <<= ((code_size-i)<<1);

            while (char_counts_[c+1] == char_counts_[c]) c++;
            parses.emplace_back(i-1, sa_[char_counts_[c]], seq[i-1]);
        } else {
            parses.emplace_back(length, sa_[a], last);
        }
    }


    uint64_t LF(const uint64_t& i) const {
        uint8_t c = at(i);
        return char_counts_[c] + rank(i, c);
    }

    uint64_t rank(uint64_t i, uint8_t c) const {
        if (char_counts_[c+1] == char_counts_[c]) return 0;
        c = alphabet_type::convert(c);
        if (i >= size_) [[unlikely]] {
            return reinterpret_cast<alphabet_type*>(p_sums_ + alphabet_type::size() * block_count_)->p_sum(c);
        }
        uint64_t s_block_i = i / SUPER_BLOCK_ELEMS;
        uint64_t res = reinterpret_cast<alphabet_type*>(p_sums_ + alphabet_type::size() * s_block_i)->p_sum(c);
        res += s_blocks_[s_block_i]->rank(c, i % SUPER_BLOCK_ELEMS);
        return res;
    }

    uint8_t operator[](size_t i) const {
        return at(i);
    }

    uint64_t interval_size(const std::string& meta_symbol, uint8_t code_size) const {
        uint8_t c = 0;

        for (uint8_t i = 0; i < code_size-1; i++) {
            if (i < meta_symbol.size()) {
                char ch = meta_symbol[i];
                c += ((ch >= 'C') + (ch >= 'G') + (ch >= 'T'));
            }
            c <<= 2;
        }
        if (code_size == meta_symbol.size()) {
            char ch = meta_symbol[code_size-1];
            c += ((ch >= 'C') + (ch >= 'G') + (ch >= 'T'));
        }

        return char_counts_[c + (1ULL<<((code_size-meta_symbol.size())*2))] - char_counts_[c];
    }

    uint64_t size() const { return size_; }
    uint64_t bytes() const { return bytes_; }

    void print() const {
        for (uint32_t i = 0; i < block_count_; i++) {
            uint64_t s = (i + 1) * SUPER_BLOCK_ELEMS;
            s = s > size_ ? size_ % SUPER_BLOCK_ELEMS : SUPER_BLOCK_ELEMS;
            std::cerr << "S block " << i << ":" << std::endl;
            reinterpret_cast<alphabet_type*>(p_sums_ * alphabet_type::size() * i)->print();
            s_blocks_[i]->print(s);
        }
    }

    private:
    super_block_type* read_super_block(std::fstream& in_file) {
        uint64_t in_bytes = 0;
        in_file.read(reinterpret_cast<char*>(&in_bytes),sizeof(uint64_t));
        uint8_t* data = static_cast<uint8_t*>( std::malloc(in_bytes + block_type::padding_bytes));

        uint64_t n_words = 0;
        in_file.read(reinterpret_cast<char*>(&n_words),sizeof(uint64_t));
        in_file.read(reinterpret_cast<char*>(data), n_words);

        uint64_t offset = 0;
        in_file.read(reinterpret_cast<char*>(&offset),sizeof(uint64_t));
        offset += n_words;

        in_file.read(reinterpret_cast<char*>(&n_words),sizeof(uint64_t));
        in_file.read(reinterpret_cast<char*>(data+offset), n_words);

        bytes_ += in_bytes;

        return reinterpret_cast<super_block_type*>(data);
    }
    };
}  // namespace bbwt
