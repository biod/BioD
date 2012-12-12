module bio.bam.md.operation;

import bio.core.base;
import bio.core.sequence;

import std.conv;
import std.traits;
import std.bitmanip;
import std.algorithm;

/// MD tag operation types
enum MdOperationType : ubyte {
    Match,
    Mismatch,
    Deletion
}

/// Single MD operation.
struct MdOperation {

    alias ReturnType!(nucleotideSequence!SliceableString) DeletionSequence;

    private {
        MdOperationType _type;
        union {
            uint _match; 
            DeletionSequence _deletion;
            Base16 _mismatch;
        }
    }

    /// Operation type
    MdOperationType type() @property const {
        return _type;
    }

    /// ditto
    void type(MdOperationType t) @property {
        _type = t;
    }

    /// Convenience methods
    bool is_deletion() @property const {
        return _type == MdOperationType.Deletion;
    }

    /// ditto
    bool is_match() @property const {
        return _type == MdOperationType.Match;
    }

    /// ditto
    bool is_mismatch() @property const {
        return _type == MdOperationType.Mismatch;
    }

    /// The number of matched bases
    ref uint match() @property {
        return _match;
    }

    /// Mismatched reference base
    Base16 mismatch() @property const {
        return _mismatch;
    }

    /// ditto
    void mismatch(Base16 base) @property {
        _mismatch = base;
    }

    /// Deleted sequence
    ref DeletionSequence deletion() @property {
        return _deletion;
    }

    static MdOperation createMatch(uint match) {
        MdOperation m;
        m._type = MdOperationType.Match;
        m._match = match;
        return m;
    }

    static MdOperation createDeletion(string deletion) {
        MdOperation m;
        m._type = MdOperationType.Deletion;
        m._deletion = nucleotideSequence(sliceableString(deletion));
        return m;
    }

    static MdOperation createMismatch(char mismatch) {
        MdOperation m;
        m._type = MdOperationType.Mismatch;
        m._mismatch = Base16(mismatch);
        return m;
    }

    static MdOperation createDeletion(DeletionSequence seq) {
        MdOperation m;
        m._type = MdOperationType.Deletion;
        m._deletion = seq;
        return m;
    }

    static MdOperation createMismatch(Base16 base) {
        MdOperation m;
        m._type = MdOperationType.Mismatch;
        m._mismatch = base;
        return m;
    }

    bool opEquals(ref MdOperation other) {

        if (type != other.type) {
            return false;
        }

        final switch (type) {
            case MdOperationType.Match:
                return _match == other._match;
            case MdOperationType.Mismatch:
                return mismatch == other.mismatch;
            case MdOperationType.Deletion:
                return equal(_deletion, other._deletion);
        }
    }

    string toString() const {
        final switch (type) {
            case MdOperationType.Match:
                return "Match(" ~ to!string(_match) ~ ")";
            case MdOperationType.Mismatch:
                return "Mismatch(" ~ to!string(_mismatch) ~ ")";
            case MdOperationType.Deletion:
                return "Deletion(" ~ to!string(_deletion) ~ ")";
        }
    }
}

/// Returns MD operation with reverse-complemented data
MdOperation reverseMdOp(MdOperation op) {
    if (op.is_deletion)
        return MdOperation.createDeletion(op.deletion.reverse);

    if (op.is_mismatch)
        return MdOperation.createMismatch(op.mismatch.complement);

    return op;
}
