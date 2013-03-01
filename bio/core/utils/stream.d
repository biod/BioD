module bio.core.utils.stream;

public import std.stream;
import core.stdc.stdio;
import core.stdc.errno;
import core.stdc.string;
import core.sys.posix.sys.select;
import std.conv;

version(Posix){
    private import core.sys.posix.unistd;
}

FileMode toFileMode(string mode) {
    FileMode result = FileMode.In;
    switch (mode) {
        case "r", "rb":
            result = FileMode.In; // 1000
            break;
        case "r+", "r+b", "rb+":
            result = FileMode.In | FileMode.Out; // 1100
        case "w", "wb":
            result = FileMode.OutNew; // 0110
            break;
        case "w+", "w+b", "wb+":
            result = FileMode.In | FileMode.OutNew; // 1110
            break;
        case "a", "ab":
            result = FileMode.Append; // 0001
            break;
        case "a+", "a+b", "ab+":
            result = FileMode.In | FileMode.Append; // 1001
            break;
        default:
            break;
    }

    return result;
}

final class File: std.stream.File {
    this(string filename, string mode="rb") {
        // Issue 8528 workaround
        auto file = fopen(toStringz(filename), toStringz(mode));
        if (file == null) {
            throw new OpenException(cast(string) ("Cannot open or create file '"
                                            ~ filename ~ "'"));
        }
        super(core.stdc.stdio.fileno(file), toFileMode(mode));
    }

    override ulong seek(long offset, SeekPos rel) {
        assertSeekable();
        auto hFile = handle();
        version (Windows) {
          int hi = cast(int)(offset>>32);
          uint low = SetFilePointer(hFile, cast(int)offset, &hi, rel);
          if ((low == INVALID_SET_FILE_POINTER) && (GetLastError() != 0))
            throw new SeekException("unable to move file pointer");
          ulong result = (cast(ulong)hi << 32) + low; 
        } else version (Posix) {
            // Phobos casts offset to int, leading to throwing an exception
            // on large files
            auto result = lseek(hFile, cast(off_t)offset, rel);
        }    
        if (result == cast(typeof(result))-1)
          throw new SeekException("unable to move file pointer");
        readEOF = false;
        return cast(ulong)result;
      }

    override size_t readBlock(void* buffer, size_t size) {
        assertReadable();
        auto hFile = handle();
        version (Windows) {
            auto dwSize = to!DWORD(size);
            ReadFile(hFile, buffer, dwSize, &dwSize, null);
            size = dwSize;
        } else version (Posix) {
            // http://developerweb.net/viewtopic.php?id=4267
            fd_set rset;
            timeval timeout;
            immutable MAX_IDLE_SECS = 1;
            while (true) {
                auto ret = core.sys.posix.unistd.read(hFile, buffer, size);
                if (ret == -1) {
                    if (errno == EINTR)
                        continue;

                    if (errno == EAGAIN || errno == EWOULDBLOCK) {
                        FD_ZERO(&rset);
                        FD_SET(hFile, &rset);
                        timeout.tv_sec = MAX_IDLE_SECS;
                        timeout.tv_usec = 0;
                        ret = select(hFile + 1, &rset, null, null, &timeout);
                        if (ret <= 0) {
                            size = 0;
                            throw new ReadException("read timeout");
                            break;
                        }
                    } else {
                        throw new ReadException(to!string(strerror(errno)));
                    }
                } else {
                    size = ret;
                    break;
                }
            }
        }
        readEOF = (size == 0);
        return size;
    }
}
