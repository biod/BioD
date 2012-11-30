module bio.core.kmer;

import bio.core.base;

///
struct KMer(uint K) {
    private ulong _id;
    ulong id() @property const {
        return _id;
    }

    this(S)(S sequence) {
        size_t i = 0;
        foreach (nuc; sequence) {
            _id <<= 2;
            ++i;
            switch (nuc) {
                case 'A':
                    break;
                case 'C':
                    _id += 1;
                    break;
                case 'G':
                    _id += 2;
                    break;
                case 'T':
                    _id += 3;
                    break;
                default:
                    _id >>= 2;
                    --i;
                    break;
            }

            if (i == K)
                break;
        }
    }

    struct KMerSequence {
        this(ulong number) {
            _n = number;
        }

        private ulong _n;
        private size_t _len = K;

        bool empty() @property const { return _len == 0; }
        void popFront() { --_len; }
        void popBack() { --_len; _n >>= 2; }

        private static Base5 code2base(int code) {
            return Base5("ACGT"[code]);
        }

        Base5 opIndex(size_t i) const {
            return code2base((_n >> (2 * (_len - i - 1))) & 3);
        }

        size_t length() @property const { return _len; }
        Base5 front() @property const { return opIndex(0); }
        Base5 back() @property const { return opIndex(_len - 1); }
        KMerSequence save() @property const { 
            KMerSequence _seq = void;
            _seq._n = _n;
            _seq._len = _len;
            return _seq;
        }
    }

    KMerSequence sequence() @property const {
        return KMerSequence(_id);
    }
}

unittest {
    import std.algorithm;
    auto kmer = KMer!10("AACGTACGTG");
    assert(equal(kmer.sequence, "AACGTACGTG"));
}
