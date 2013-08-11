module bio.core.utils.file;

import std.stdio;
import std.system;
import std.bitmanip;
import std.array;

bool isSeekable(ref File f) @property {
    ulong pos = void;
    bool result = false;
    try {
        pos = f.tell;
        result = fseek(f.getFP(), 0, 0) != -1;
    } catch (Exception e) {
        return result;
    }
    f.seek(pos);
    return result;
}

void seekCur(ref File f, ulong offset) {
    f.seek(offset, SEEK_CUR);
}

auto peekLE(T, R)(R range) {
    return range.peek!(T, Endian.littleEndian)();
}

auto peekBE(T, R)(R range) {
    return range.peek!(T, Endian.bigEndian)();
}
