/**
 * Contains the obsolete pattern matching functions from Phobos'
 * `std.string`.
 */
module undead.string;

import std.traits;

/***********************************************
 * See if character c is in the pattern.
 * Patterns:
 *
 *  A $(I pattern) is an array of characters much like a $(I character
 *  class) in regular expressions. A sequence of characters
 *  can be given, such as "abcde". The '-' can represent a range
 *  of characters, as "a-e" represents the same pattern as "abcde".
 *  "a-fA-F0-9" represents all the hex characters.
 *  If the first character of a pattern is '^', then the pattern
 *  is negated, i.e. "^0-9" means any character except a digit.
 *  The functions inPattern, $(B countchars), $(B removeschars),
 *  and $(B squeeze) use patterns.
 *
 * Note: In the future, the pattern syntax may be improved
 *  to be more like regular expression character classes.
 */
bool inPattern(S)(dchar c, in S pattern) @safe pure @nogc
if (isSomeString!S)
{
    bool result = false;
    int range = 0;
    dchar lastc;

    foreach (size_t i, dchar p; pattern)
    {
        if (p == '^' && i == 0)
        {
            result = true;
            if (i + 1 == pattern.length)
                return (c == p);    // or should this be an error?
        }
        else if (range)
        {
            range = 0;
            if (lastc <= c && c <= p || c == p)
                return !result;
        }
        else if (p == '-' && i > result && i + 1 < pattern.length)
        {
            range = 1;
            continue;
        }
        else if (c == p)
            return !result;
        lastc = p;
    }
    return result;
}


@safe pure @nogc unittest
{
    assertCTFEable!(
    {
    assert(inPattern('x', "x") == 1);
    assert(inPattern('x', "y") == 0);
    assert(inPattern('x', string.init) == 0);
    assert(inPattern('x', "^y") == 1);
    assert(inPattern('x', "yxxy") == 1);
    assert(inPattern('x', "^yxxy") == 0);
    assert(inPattern('x', "^abcd") == 1);
    assert(inPattern('^', "^^") == 0);
    assert(inPattern('^', "^") == 1);
    assert(inPattern('^', "a^") == 1);
    assert(inPattern('x', "a-z") == 1);
    assert(inPattern('x', "A-Z") == 0);
    assert(inPattern('x', "^a-z") == 0);
    assert(inPattern('x', "^A-Z") == 1);
    assert(inPattern('-', "a-") == 1);
    assert(inPattern('-', "^A-") == 0);
    assert(inPattern('a', "z-a") == 1);
    assert(inPattern('z', "z-a") == 1);
    assert(inPattern('x', "z-a") == 0);
    });
}


/**
 * See if character c is in the intersection of the patterns.
 */
bool inPattern(S)(dchar c, S[] patterns) @safe pure @nogc
if (isSomeString!S)
{
    foreach (string pattern; patterns)
    {
        if (!inPattern(c, pattern))
        {
            return false;
        }
    }
    return true;
}


/**
 * Count characters in s that match pattern.
 */
size_t countchars(S, S1)(S s, in S1 pattern) @safe pure @nogc
if (isSomeString!S && isSomeString!S1)
{
    size_t count;
    foreach (dchar c; s)
    {
        count += inPattern(c, pattern);
    }
    return count;
}

@safe pure @nogc unittest
{
    assertCTFEable!(
    {
    assert(countchars("abc", "a-c") == 3);
    assert(countchars("hello world", "or") == 3);
    });
}


/**
 * Return string that is s with all characters removed that match pattern.
 */
S removechars(S)(S s, in S pattern) @safe pure
if (isSomeString!S)
{
    import std.utf : encode;

    Unqual!(typeof(s[0]))[] r;
    bool changed = false;

    foreach (size_t i, dchar c; s)
    {
        if (inPattern(c, pattern))
        {
            if (!changed)
            {
                changed = true;
                r = s[0 .. i].dup;
            }
            continue;
        }
        if (changed)
        {
            encode(r, c);
        }
    }
    if (changed)
        return r;
    else
        return s;
}

@safe pure unittest
{
    assertCTFEable!(
    {
    assert(removechars("abc", "a-c").length == 0);
    assert(removechars("hello world", "or") == "hell wld");
    assert(removechars("hello world", "d") == "hello worl");
    assert(removechars("hah", "h") == "a");
    });
}

@safe pure unittest
{
    assert(removechars("abc", "x") == "abc");
}


/***************************************************
 * Return string where sequences of a character in s[] from pattern[]
 * are replaced with a single instance of that character.
 * If pattern is null, it defaults to all characters.
 */
S squeeze(S)(S s, in S pattern = null)
{
    import std.utf : encode, stride;

    Unqual!(typeof(s[0]))[] r;
    dchar lastc;
    size_t lasti;
    int run;
    bool changed;

    foreach (size_t i, dchar c; s)
    {
        if (run && lastc == c)
        {
            changed = true;
        }
        else if (pattern is null || inPattern(c, pattern))
        {
            run = 1;
            if (changed)
            {
                if (r is null)
                    r = s[0 .. lasti].dup;
                encode(r, c);
            }
            else
                lasti = i + stride(s, i);
            lastc = c;
        }
        else
        {
            run = 0;
            if (changed)
            {
                if (r is null)
                    r = s[0 .. lasti].dup;
                encode(r, c);
            }
        }
    }
    return changed ? ((r is null) ? s[0 .. lasti] : cast(S) r) : s;
}

@system pure unittest
{
    assertCTFEable!(
    {
    string s;

    assert(squeeze("hello") == "helo");

    s = "abcd";
    assert(squeeze(s) is s);
    s = "xyzz";
    assert(squeeze(s).ptr == s.ptr); // should just be a slice

    assert(squeeze("hello goodbyee", "oe") == "hello godbye");
    });
}

/***************************************************************
 Finds the position $(D_PARAM pos) of the first character in $(D_PARAM
 s) that does not match $(D_PARAM pattern) (in the terminology used by
 $(REF inPattern, std,string)). Updates $(D_PARAM s =
 s[pos..$]). Returns the slice from the beginning of the original
 (before update) string up to, and excluding, $(D_PARAM pos).

The $(D_PARAM munch) function is mostly convenient for skipping
certain category of characters (e.g. whitespace) when parsing
strings. (In such cases, the return value is not used.)
 */
S1 munch(S1, S2)(ref S1 s, S2 pattern) @safe pure @nogc
{
    size_t j = s.length;
    foreach (i, dchar c; s)
    {
        if (!inPattern(c, pattern))
        {
            j = i;
            break;
        }
    }
    scope(exit) s = s[j .. $];
    return s[0 .. j];
}

///
@safe pure @nogc unittest
{
    string s = "123abc";
    string t = munch(s, "0123456789");
    assert(t == "123" && s == "abc");
    t = munch(s, "0123456789");
    assert(t == "" && s == "abc");
}

@safe pure @nogc unittest
{
    string s = "123€abc";
    string t = munch(s, "0123456789");
    assert(t == "123" && s == "€abc");
    t = munch(s, "0123456789");
    assert(t == "" && s == "€abc");
    t = munch(s, "£$€¥");
    assert(t == "€" && s == "abc");
}

// helper function for unit tests
private @property void assertCTFEable(alias dg)()
{
    static assert({ cast(void) dg(); return true; }());
    cast(void) dg();
}