// Build: g++ -std=gnu++17 -O2 -pipe -static -s main.cpp -o solve
// Run:   ./solve < input.json
// Output: prints secret (P(0)) and indices of wrong shares

#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>

using namespace std;
using boost::multiprecision::cpp_int;

// Simple big rational with gcd reduction
struct Fraction {
    cpp_int num; // numerator
    cpp_int den; // denominator > 0
    Fraction(): num(0), den(1) {}
    Fraction(const cpp_int& n): num(n), den(1) {}
    Fraction(const cpp_int& n, const cpp_int& d) { set(n,d); }

    static cpp_int igcd(cpp_int a, cpp_int b){
        if (a < 0) a = -a;
        if (b < 0) b = -b;
        while (b != 0) { cpp_int t = a % b; a = b; b = t; }
        return a;
    }

    void normalize() {
        if (den < 0) { den = -den; num = -num; }
        if (num == 0) { den = 1; return; }
        cpp_int g = igcd(num, den);
        if (g != 0) { num /= g; den /= g; }
    }
    void set(const cpp_int& n, const cpp_int& d) {
        num = n; den = d; normalize();
    }

    static Fraction add(const Fraction& a, const Fraction& b) {
        Fraction r;
        r.num = a.num * b.den + b.num * a.den;
        r.den = a.den * b.den;
        r.normalize();
        return r;
    }
    static Fraction sub(const Fraction& a, const Fraction& b) {
        Fraction r;
        r.num = a.num * b.den - b.num * a.den;
        r.den = a.den * b.den;
        r.normalize();
        return r;
    }
    static Fraction mul(const Fraction& a, const Fraction& b) {
        Fraction r;
        r.num = a.num * b.num;
        r.den = a.den * b.den;
        r.normalize();
        return r;
    }
    static Fraction div(const Fraction& a, const Fraction& b) {
        Fraction r;
        r.num = a.num * b.den;
        r.den = a.den * b.num;
        r.normalize();
        return r;
    }

    bool isInteger() const { return den == 1; }
    bool operator==(const Fraction& o) const { return num == o.num && den == o.den; }
};

struct Share {
    int x;          // x-coordinate (index)
    int base;       // base of value string
    string sval;    // value string (y in given base)
    cpp_int y;      // parsed y (big int)
};

int digitVal(char c){
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'z') return 10 + (c - 'a');
    if (c >= 'A' && c <= 'Z') return 10 + (c - 'A');
    return -1;
}

// parse base-b string to cpp_int (b up to 36 for generality)
cpp_int parseBase(const string& s, int base){
    cpp_int val = 0;
    for (char c : s) {
        if (c == '_' || c == ' ') continue;
        int d = digitVal(c);
        if (d < 0 || d >= base) {
            throw runtime_error("Invalid digit for base");
        }
        val = val * base + d;
    }
    return val;
}

// Lagrange basis lambda_i(0) = prod_{j != i} (-x_j)/(x_i - x_j)
Fraction lagrange_lambda_at_zero(const vector<int>& xs, int i){
    Fraction lam( cpp_int(1) );
    for (int j = 0; j < (int)xs.size(); ++j){
        if (j == i) continue;
        cpp_int num = -cpp_int(xs[j]);
        cpp_int den = cpp_int(xs[i]) - cpp_int(xs[j]);
        lam = Fraction::mul(lam, Fraction(num, den));
    }
    return lam;
}

// P(0) from subset using Lagrange form
Fraction interpolate_P0(const vector<int>& xs, const vector<cpp_int>& ys){
    int k = xs.size();
    Fraction sum( cpp_int(0) );
    for (int i = 0; i < k; ++i){
        Fraction lam = lagrange_lambda_at_zero(xs, i);
        Fraction term = Fraction::mul(Fraction(ys[i]), lam);
        sum = Fraction::add(sum, term);
    }
    return sum;
}

// Evaluate P at query xq using Lagrange basis
Fraction interpolate_eval(const vector<int>& xs, const vector<cpp_int>& ys, int xq){
    int k = xs.size();
    // If xq equals one of xs, return exact y to avoid 0/0 numerics
    for (int i = 0; i < k; ++i){
        if (xq == xs[i]) return Fraction(ys[i]);
    }
    Fraction sum( cpp_int(0) );
    for (int i = 0; i < k; ++i){
        // l_i(xq) = prod_{j != i} (xq - x_j)/(x_i - x_j)
        Fraction li( cpp_int(1) );
        for (int j = 0; j < k; ++j){
            if (i == j) continue;
            cpp_int num = cpp_int(xq) - cpp_int(xs[j]);
            cpp_int den = cpp_int(xs[i]) - cpp_int(xs[j]);
            li = Fraction::mul(li, Fraction(num, den));
        }
        Fraction term = Fraction::mul(Fraction(ys[i]), li);
        sum = Fraction::add(sum, term);
    }
    return sum;
}

// Simple hand-rolled minimal JSON reader tailored to given format
// Assumes well-formed input like the prompt; robust parsing is not the goal here.
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Read entire stdin
    string json((istreambuf_iterator<char>(cin)), istreambuf_iterator<char>());

    // Extract n and k
    auto findNumberAfter = [&](const string& key)->long long{
        size_t p = json.find("\"" + key + "\"");
        if (p == string::npos) throw runtime_error("Missing key");
        p = json.find(':', p);
        if (p == string::npos) throw runtime_error("Malformed");
        ++p;
        while (p < json.size() && isspace((unsigned char)json[p])) ++p;
        // read number
        long long val = 0;
        bool neg = false;
        if (json[p] == '-') { neg = true; ++p; }
        if (!isdigit((unsigned char)json[p])) throw runtime_error("Expected digit");
        while (p < json.size() && isdigit((unsigned char)json[p])) {
            val = val * 10 + (json[p]-'0');
            ++p;
        }
        return neg ? -val : val;
    };

    int n = (int)findNumberAfter("n");
    int k = (int)findNumberAfter("k");

    vector<Share> shares;
    shares.reserve(n);

    // For indices 1..n, try to find base and value
    for (int idx = 1; idx <= n; ++idx){
        string key = "\"" + to_string(idx) + "\"";
        size_t p = json.find(key);
        if (p == string::npos) {
            // This index absent; skip
            continue;
        }
        // find base
        size_t pb = json.find("\"base\"", p);
        if (pb == string::npos) throw runtime_error("Missing base");
        pb = json.find(':', pb);
        if (pb == string::npos) throw runtime_error("Malformed base");
        ++pb;
        while (pb < json.size() && (isspace((unsigned char)json[pb]) || json[pb]=='\"')) ++pb;
        // read base number up to quote or comma
        string baseTok;
        while (pb < json.size() && isalnum((unsigned char)json[pb])) {
            baseTok.push_back(json[pb++]);
        }
        int base = stoi(baseTok);

        // find value
        size_t pv = json.find("\"value\"", p);
        if (pv == string::npos) throw runtime_error("Missing value");
        pv = json.find(':', pv);
        if (pv == string::npos) throw runtime_error("Malformed value");
        pv = json.find('\"', pv);
        if (pv == string::npos) throw runtime_error("Missing value quote");
        ++pv;
        string valstr;
        while (pv < json.size() && json[pv] != '\"') {
            valstr.push_back(json[pv++]);
        }

        Share sh;
        sh.x = idx;
        sh.base = base;
        sh.sval = valstr;
        sh.y = parseBase(valstr, base);
        shares.push_back(std::move(sh));
    }

    if ((int)shares.size() < k) {
        cerr << "Not enough shares to reconstruct\n";
        return 1;
    }

    // Prepare data containers
    // For combinatorics, collect available indices (by position in shares vector)
    vector<int> idxs(shares.size());
    iota(idxs.begin(), idxs.end(), 0);

    // Choose best combination by inliers
    cpp_int bestSecret = 0;
    int bestInliers = -1;
    vector<int> bestSubset; // positions
    vector<int> bestInlierMask;

    // Generate all k-combinations using bitmask lexicographic
    vector<int> select;
    function<void(int,int)> dfs = [&](int start, int need){
        static vector<int> cur;
        if ((int)cur.size() + (int)idxs.size() - start < need) return;
        if (need == 0) {
            // Build xs, ys from cur
            vector<int> xs;
            vector<cpp_int> ys;
            xs.reserve(k);
            ys.reserve(k);
            for (int pos : cur) {
                xs.push_back(shares[pos].x);
                ys.push_back(shares[pos].y);
            }
            // Compute P(0)
            Fraction p0 = interpolate_P0(xs, ys);
            if (!p0.isInteger()) {
                // skip non-integer secret given problem expectation
                return;
            }
            cpp_int secret = p0.num; // integer

            // Count inliers: evaluate at each share.x and compare to y
            int inliers = 0;
            vector<int> inlierMask(shares.size(), 0);
            for (int t = 0; t < (int)shares.size(); ++t){
                Fraction val = interpolate_eval(xs, ys, shares[t].x);
                if (val.den == 0) continue; // shouldn't happen
                if (val.den != 1) continue; // must be integer y
                if (val.num == shares[t].y) {
                    inliers++;
                    inlierMask[t] = 1;
                }
            }
            if (inliers > bestInliers) {
                bestInliers = inliers;
                bestSecret = secret;
                bestSubset = cur;
                bestInlierMask = inlierMask;
            }
            return;
        }
        for (int i = start; i < (int)idxs.size(); ++i){
            cur.push_back(idxs[i]);
            dfs(i+1, need-1);
            cur.pop_back();
        }
    };
    dfs(0, k);

    // Output
    cout << "secret=" << bestSecret << "\n";
    cout << "wrong_share_indices=[";
    bool first = true;
    for (int i = 0; i < (int)shares.size(); ++i){
        if (!bestInlierMask.empty() && bestInlierMask[i] == 0) {
            if (!first) cout << ",";
            cout << shares[i].x;
            first = false;
        }
    }
    cout << "]\n";

    return 0;
}

