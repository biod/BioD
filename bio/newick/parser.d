module bio.newick.parser;

import bio.newick.treenode;

import std.exception;
import std.format;
import std.ascii;
import std.array;

class NewickParseException : Exception {
    this(string msg) {
        super(msg);
    }
}

NewickNode* parse(string s) {
    string str = s;
    return parseTree(str);
}

private {

    void skipWhiteSpace(ref string s) {
        while (true) {
            if (s.empty)
                throw new NewickParseException("String is empty");

            if (isWhite(s.front))
                s = s[1 .. $];
            else
                break;
        }
    }

    NewickNode* parseTree(ref string s) {

        skipWhiteSpace(s);
        auto ch = s.front;

        if (ch == ';')
            return new NewickNode;

        auto result = parseSubtree(s);
        enforce(s.front == ';', "String does not end with ';'");
        return result;
    }

    NewickNode* parseSubtree(ref string s) {
        skipWhiteSpace(s);

        auto ch = s.front;

        switch (ch) {
            case '(':
                return parseInternalNode(s);

            case ',':
            case ')':
                return new NewickNode;

            default:
                enforce(isAlphaNum(ch), "Node name contains invalid character");
                return parseLeafNode(s);
        }
    }

    NewickNode* parseLeafNode(ref string s, NewickNode* nd=null) {
        skipWhiteSpace(s);

        auto ch = s.front;

        auto node = nd is null ? (new NewickNode) : nd;

        size_t i = 0;
        while (isAlphaNum(s[i])) {
            ++i;
        }

        if (i > 0)
            node.name = s[0 .. i];
        
        s = s[i .. $];

        ch = s.front;

        if (ch == ':')
            formattedRead(s, ":%s", &node.distance_to_parent);

        return node;
    }

    NewickNode* parseInternalNode(ref string s) {
        skipWhiteSpace(s);

        auto ch = s.front;
        assert(ch == '(');
        s.popFront();

        auto node = new NewickNode;

        auto children = Appender!(NewickNode*[])();
            
        while (true) {
            auto child = parseSubtree(s);
            children.put(child);

            ch = s.front;

            if (ch == ',') {
                s.popFront();
                continue;
            } else if (ch == ')') {
                s.popFront();
                break;
            } else {
                throw new NewickParseException("Parse error");
            }
        }

        node.children = children.data;
        parseLeafNode(s, node);

        return node;
    }
}
