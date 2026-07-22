#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kAvogadroScale = 0.602; // Converts (g/mol)/(g/cm^3) to A^3.

struct Settings {
    int n1 = 128, m1 = 900; // Bifunctional network strands.
    int n2 = 32,  m2 = 225; // Multifunctional cross-linkers.
    int n3 = 32,  m3 = 425; // Neutral linear oil fillers.
    int n4 = 5,   m4 = 6;   // Five-bead star-like moderators.
    int crosslinker_functionality = 8;
    std::string crosslink_distribution = "random";
    std::uint32_t crosslink_seed = 20260722u;
    double bead_mass = 74.0;
    double density = 0.1;
    double bond_length = 2.801;
    double spacing = 7.0;
    std::uint32_t seed = 5489u;
    std::string output = "data.V22_PDMS_N32_10wt";
    bool output_explicit = false;
    bool filler_length_explicit = false;
    double filler_weight_percent = -1.0;
};

struct Atom { int id, molecule, type; double charge, x, y, z; };
struct Bond { int id, type, a, b; };
struct Angle { int id, type, a, b, c; };
struct Dihedral { int id, type, a, b, c, d; };

struct System {
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::vector<Angle> angles;
    std::vector<Dihedral> dihedrals;
};

int parse_int(const std::string& value, const std::string& option) {
    try { size_t used = 0; int result = std::stoi(value, &used); if (used != value.size()) throw std::invalid_argument(""); return result; }
    catch (...) { throw std::runtime_error("Invalid integer for " + option + ": " + value); }
}

double parse_double(const std::string& value, const std::string& option) {
    try { size_t used = 0; double result = std::stod(value, &used); if (used != value.size()) throw std::invalid_argument(""); return result; }
    catch (...) { throw std::runtime_error("Invalid number for " + option + ": " + value); }
}

void print_help(const char* program) {
    std::cout
        << "Usage: " << program << " [options]\n\n"
        << "V22 LAMMPS data-file generator. N is beads per molecule; M is molecule count.\n\n"
        << "  --n1 N --m1 M    bifunctional network strands (defaults: 128, 900)\n"
        << "  --n2 N           cross-linker length (default: 32); M2 is stoichiometric\n"
        << "  --n3 N --m3 M    neutral linear oil fillers (defaults: 32, 425)\n"
        << "  --n4 N --m4 M    star-moderator size/count (N4 must be 5; defaults: 5, 6)\n"
        << "  --filler-length N  PDMS filler chain length (alias for --n3)\n"
        << "  --filler-wt X      PDMS weight percent; computes M3 from fixed V22 base\n"
        << "  --functionality F cross-linker reactive sites, 3-16 (default: 8)\n"
        << "  --crosslink-distribution MODE\n"
        << "                     reactive-site placement: regular or random (default: random)\n"
        << "  --crosslink-seed N random-site seed (default: 20260722)\n"
        << "  --mass X          bead mass in g/mol (default: 74)\n"
        << "  --density X       target density in g/cm^3 (default: 0.1)\n"
        << "  --bond-length X   initial bond length in angstrom (default: 2.801)\n"
        << "  --spacing X       spacing between placed molecules (default: 7.0)\n"
        << "  --seed N          reproducible random seed (default: 5489)\n"
        << "  --output FILE     override the automatically generated output filename\n"
        << "  --help             show this help\n";
}

Settings parse_args(int argc, char** argv) {
    Settings s;
    for (int i = 1; i < argc; ++i) {
        const std::string option = argv[i];
        if (option == "--help") { print_help(argv[0]); std::exit(0); }
        if (i + 1 >= argc) throw std::runtime_error("Missing value after " + option);
        const std::string value = argv[++i];
        if      (option == "--n1") s.n1 = parse_int(value, option);
        else if (option == "--m1") s.m1 = parse_int(value, option);
        else if (option == "--n2") s.n2 = parse_int(value, option);
        else if (option == "--m2")
            throw std::runtime_error("M2 is determined by M2=2*M1/functionality; do not supply --m2");
        else if (option == "--n3") s.n3 = parse_int(value, option);
        else if (option == "--m3") s.m3 = parse_int(value, option);
        else if (option == "--n4") s.n4 = parse_int(value, option);
        else if (option == "--m4") s.m4 = parse_int(value, option);
        else if (option == "--functionality") s.crosslinker_functionality = parse_int(value, option);
        else if (option == "--crosslink-distribution") s.crosslink_distribution = value;
        else if (option == "--crosslink-seed") s.crosslink_seed = static_cast<std::uint32_t>(parse_int(value, option));
        else if (option == "--mass") s.bead_mass = parse_double(value, option);
        else if (option == "--density") s.density = parse_double(value, option);
        else if (option == "--bond-length") s.bond_length = parse_double(value, option);
        else if (option == "--spacing") s.spacing = parse_double(value, option);
        else if (option == "--seed") s.seed = static_cast<std::uint32_t>(parse_int(value, option));
        else if (option == "--output") { s.output = value; s.output_explicit = true; }
        else if (option == "--filler-length") { s.n3 = parse_int(value, option); s.filler_length_explicit = true; }
        else if (option == "--filler-wt") s.filler_weight_percent = parse_double(value, option);
        else throw std::runtime_error("Unknown option: " + option);
    }
    return s;
}

std::string filename_number(double value) {
    std::ostringstream text;
    text << std::fixed << std::setprecision(4) << value;
    std::string result = text.str();
    while (!result.empty() && result.back() == '0') result.pop_back();
    if (!result.empty() && result.back() == '.') result.pop_back();
    std::replace(result.begin(), result.end(), '.', 'p');
    return result;
}

void apply_filler_composition(Settings& s) {
    if (s.filler_length_explicit && s.filler_weight_percent < 0.0)
        s.filler_weight_percent = 10.0;
    if (s.filler_weight_percent < 0.0) return;
    if (s.filler_weight_percent >= 100.0)
        throw std::runtime_error("--filler-wt must be at least 0 and less than 100");
    if (s.n3 <= 0 && s.filler_weight_percent > 0.0)
        throw std::runtime_error("A positive --filler-wt requires a positive --filler-length");

    // Components 1, 2, and 4 form the fixed V22 base. With equal bead masses,
    // w = filler_beads / (base_beads + filler_beads).
    const long long base_beads = 1LL*s.n1*s.m1 + 1LL*s.n2*s.m2 + 1LL*s.n4*s.m4;
    if (base_beads <= 0) throw std::runtime_error("The fixed V22 base must contain at least one bead");
    if (s.filler_weight_percent == 0.0) {
        s.m3 = 0;
    } else {
        const double fraction = s.filler_weight_percent / 100.0;
        const double desired_filler_beads = base_beads * fraction / (1.0 - fraction);
        s.m3 = std::max(1, static_cast<int>(std::lround(desired_filler_beads / s.n3)));
    }
    if (!s.output_explicit)
        s.output = "data.V22_PDMS_N" + std::to_string(s.n3) + "_" + filename_number(s.filler_weight_percent) + "wt";
}

void apply_crosslinker_stoichiometry(Settings& s) {
    if (s.crosslinker_functionality < 3)
        throw std::runtime_error("Cross-linker functionality must be at least 3 to form a network");
    const long long reactive_ends = 2LL * s.m1;
    if (reactive_ends % s.crosslinker_functionality != 0) {
        std::ostringstream message;
        message << "Exact stoichiometry is impossible: 2*M1=" << reactive_ends
                << " is not divisible by functionality=" << s.crosslinker_functionality;
        throw std::runtime_error(message.str());
    }
    const long long molecule_count = reactive_ends / s.crosslinker_functionality;
    if (molecule_count > 100000000LL)
        throw std::runtime_error("Stoichiometric M2 is too large");
    s.m2 = static_cast<int>(molecule_count);
}

void report_composition(const Settings& s) {
    const int n[] = {s.n1, s.n2, s.n3, s.n4};
    const int m[] = {s.m1, s.m2, s.m3, s.m4};
    long long total_beads = 0;
    long long total_molecules = 0;
    for (size_t i = 0; i < 4; ++i) total_beads += 1LL * n[i] * m[i];
    for (int count : m) total_molecules += count;
    std::cerr << "Composition (equal bead masses):\n";
    for (size_t i = 0; i < 4; ++i) {
        const long long component_beads = 1LL * n[i] * m[i];
        const double realized = total_beads == 0 ? 0.0 : 100.0 * component_beads / total_beads;
        const double mole_fraction = total_molecules == 0 ? 0.0 : 100.0 * m[i] / total_molecules;
        std::cerr << "  component " << i + 1 << ": N=" << n[i] << ", M=" << m[i]
                  << ", mole%=" << std::fixed << std::setprecision(4) << mole_fraction
                  << ", realized wt%=" << realized;
        if (i == 2 && s.filler_weight_percent >= 0.0)
            std::cerr << ", requested filler wt%=" << s.filler_weight_percent;
        std::cerr << '\n';
    }
    std::cerr << "  total beads=" << total_beads << '\n';
}

void validate(const Settings& s) {
    const int values[] = {s.n1, s.n2, s.n3, s.n4, s.m1, s.m2, s.m3, s.m4};
    if (std::any_of(std::begin(values), std::end(values), [](int x) { return x < 0; }))
        throw std::runtime_error("N and M values cannot be negative");
    if (s.n4 != 5 && s.m4 != 0)
        throw std::runtime_error("The implemented moderator is a five-bead star; use --n4 5");
    if (s.crosslinker_functionality < 3 || s.crosslinker_functionality > 16 || s.crosslinker_functionality > s.n2)
        throw std::runtime_error("Cross-linker functionality must be between 3 and min(16, N2)");
    if (s.crosslink_distribution != "regular" && s.crosslink_distribution != "random")
        throw std::runtime_error("--crosslink-distribution must be regular or random");
    if (s.crosslink_distribution == "regular" && s.crosslinker_functionality > (s.n2 + 1) / 2)
        throw std::runtime_error("Regular placement at 1,3,5,... requires functionality <= ceil(N2/2)");
    if (s.bead_mass <= 0 || s.density <= 0 || s.bond_length <= 0 || s.spacing <= 0)
        throw std::runtime_error("Mass, density, bond length, and spacing must be positive");
    const long long total_beads = 1LL*s.n1*s.m1 + 1LL*s.n2*s.m2 + 1LL*s.n3*s.m3 + 1LL*s.n4*s.m4;
    if (total_beads <= 0) throw std::runtime_error("The formulation must contain at least one bead");
    if (s.output.empty()) throw std::runtime_error("Output filename cannot be empty");
}

std::set<int> regular_sites(int functionality) {
    std::set<int> sites;
    for (int i = 0; i < functionality; ++i) sites.insert(1 + 2*i);
    return sites;
}

std::set<int> random_sites(int bead_count, int functionality, std::uint32_t seed) {
    std::vector<int> candidates(static_cast<size_t>(bead_count));
    for (int i = 0; i < bead_count; ++i) candidates[static_cast<size_t>(i)] = i + 1;
    std::mt19937 rng(seed);
    std::shuffle(candidates.begin(), candidates.end(), rng);
    return std::set<int>(candidates.begin(), candidates.begin() + functionality);
}

std::set<int> crosslinker_sites(const Settings& s) {
    if (s.crosslink_distribution == "regular")
        return regular_sites(s.crosslinker_functionality);
    return random_sites(s.n2, s.crosslinker_functionality, s.crosslink_seed);
}

void add_linear_topology(System& sys, int first_atom, int bead_count) {
    for (int i = 0; i + 1 < bead_count; ++i)
        sys.bonds.push_back({static_cast<int>(sys.bonds.size()) + 1, 1, first_atom + i, first_atom + i + 1});
    for (int i = 0; i + 2 < bead_count; ++i)
        sys.angles.push_back({static_cast<int>(sys.angles.size()) + 1, 1, first_atom + i, first_atom + i + 1, first_atom + i + 2});
    for (int i = 0; i + 3 < bead_count; ++i)
        sys.dihedrals.push_back({static_cast<int>(sys.dihedrals.size()) + 1, 1, first_atom + i, first_atom + i + 1, first_atom + i + 2, first_atom + i + 3});
}

void add_linear_component(System& sys, int bead_count, int molecule_count, int component,
                          double box, double bond_length, double spacing,
                          const std::set<int>& reactive_sites = {}) {
    if (bead_count == 0 && molecule_count != 0) throw std::runtime_error("A nonzero molecule count requires N > 0");
    if (molecule_count == 0) return;

    const double span = std::max(0, bead_count - 1) * bond_length;
    const int nx = std::max(1, static_cast<int>(box / (span + spacing)));
    const int ny = std::max(1, static_cast<int>(box / spacing));
    const long long capacity = 1LL * nx * ny * ny;
    if (capacity < molecule_count) {
        std::ostringstream msg;
        msg << "Placement grid capacity (" << capacity << ") is smaller than M" << component
            << " (" << molecule_count << "). Reduce spacing or molecule length.";
        throw std::runtime_error(msg.str());
    }

    for (int molecule_index = 0; molecule_index < molecule_count; ++molecule_index) {
        const int ix = molecule_index % nx;
        const int iy = (molecule_index / nx) % ny;
        const int iz = molecule_index / (nx * ny);
        const bool reverse = component != 1;
        double x = reverse ? box / 2 - spacing - ix * (span + spacing)
                           : -box / 2 + spacing + ix * (span + spacing);
        const double y = (reverse ? box / 2 - spacing : -box / 2 + spacing) + (reverse ? -1 : 1) * iy * spacing;
        double z = 0.0;
        if (component == 1) z = -box / 2 + spacing + iz * spacing;
        else if (component == 2) z = box / 2 - spacing - iz * spacing;
        else z = -spacing - iz * spacing;

        const int molecule_id = static_cast<int>(sys.atoms.empty() ? 1 : sys.atoms.back().molecule + 1);
        const int first_atom = static_cast<int>(sys.atoms.size()) + 1;
        for (int bead = 1; bead <= bead_count; ++bead) {
            int atom_type = 1;
            if (component == 1 && (bead == 1 || bead == bead_count)) atom_type = 2;
            if (component == 2 && reactive_sites.count(bead)) atom_type = 3;
            sys.atoms.push_back({static_cast<int>(sys.atoms.size()) + 1, molecule_id, atom_type, 0.0, x, y, z});
            x += reverse ? -bond_length : bond_length;
        }
        add_linear_topology(sys, first_atom, bead_count);
    }
}

void add_star_moderators(System& sys, const Settings& s, double box, std::mt19937& rng) {
    std::uniform_real_distribution<double> xy(-box / 2, box / 2);
    std::uniform_real_distribution<double> zdist(0.02 * box, 0.12 * box); // Preserves V22's original slab.
    for (int i = 0; i < s.m4; ++i) {
        const int molecule_id = static_cast<int>(sys.atoms.empty() ? 1 : sys.atoms.back().molecule + 1);
        const int center = static_cast<int>(sys.atoms.size()) + 1;
        const double x = xy(rng), y = xy(rng), z = zdist(rng), b = s.bond_length;
        // Preserve the legacy force-field mapping: ordinary center, type-2 arms.
        sys.atoms.push_back({center,     molecule_id, 1, 0.0, x,     y,     z});
        sys.atoms.push_back({center + 1, molecule_id, 2, 0.0, x + b, y,     z});
        sys.atoms.push_back({center + 2, molecule_id, 2, 0.0, x - b, y,     z});
        sys.atoms.push_back({center + 3, molecule_id, 2, 0.0, x,     y + b, z});
        sys.atoms.push_back({center + 4, molecule_id, 2, 0.0, x,     y - b, z});
        for (int arm = 1; arm <= 4; ++arm)
            sys.bonds.push_back({static_cast<int>(sys.bonds.size()) + 1, 2, center, center + arm});
        for (int a = 1; a <= 4; ++a)
            for (int c = a + 1; c <= 4; ++c)
                sys.angles.push_back({static_cast<int>(sys.angles.size()) + 1, 1, center + a, center, center + c});
    }
}

System build_system(const Settings& s, double box) {
    System sys;
    const long long expected_atoms = 1LL * s.n1*s.m1 + 1LL * s.n2*s.m2 + 1LL * s.n3*s.m3 + 1LL * s.n4*s.m4;
    if (expected_atoms > 100000000LL) throw std::runtime_error("Requested system is too large");
    sys.atoms.reserve(static_cast<size_t>(expected_atoms));
    add_linear_component(sys, s.n1, s.m1, 1, box, s.bond_length, s.spacing);
    const std::set<int> reactive_sites = crosslinker_sites(s);
    std::cerr << "Cross-linker type-3 sites (" << s.crosslink_distribution << "):";
    for (int site : reactive_sites) std::cerr << ' ' << site;
    std::cerr << '\n';
    add_linear_component(sys, s.n2, s.m2, 2, box, s.bond_length, s.spacing,
                         reactive_sites);
    add_linear_component(sys, s.n3, s.m3, 3, box, s.bond_length, s.spacing);
    std::mt19937 rng(s.seed);
    add_star_moderators(sys, s, box, rng);
    if (static_cast<long long>(sys.atoms.size()) != expected_atoms)
        throw std::logic_error("Internal error: generated atom count does not match requested composition");
    return sys;
}

void write_data(const Settings& s, const System& sys, double box) {
    std::ofstream out(s.output);
    if (!out) throw std::runtime_error("Cannot open output file: " + s.output);
    out << "LAMMPS data file for the V22 four-component formulation\n\n"
        << sys.atoms.size() << " atoms\n4 atom types\n"
        << sys.bonds.size() << " bonds\n2 bond types\n"
        << sys.angles.size() << " angles\n1 angle types\n"
        << sys.dihedrals.size() << " dihedrals\n1 dihedral types\n\n"
        << std::fixed << std::setprecision(6)
        << -box/2 << ' ' << box/2 << " xlo xhi\n"
        << -box/2 << ' ' << box/2 << " ylo yhi\n"
        << -box/2 << ' ' << box/2 << " zlo zhi\n\nMasses\n\n";
    for (int type = 1; type <= 4; ++type) out << type << ' ' << s.bead_mass << '\n';
    out << "\nAtoms # full\n\n" << std::setprecision(3);
    for (const auto& a : sys.atoms)
        out << a.id << ' ' << a.molecule << ' ' << a.type << ' ' << a.charge << ' '
            << a.x << ' ' << a.y << ' ' << a.z << '\n';
    if (!sys.bonds.empty()) {
        out << "\nBonds\n\n";
        for (const auto& b : sys.bonds) out << b.id << ' ' << b.type << ' ' << b.a << ' ' << b.b << '\n';
    }
    if (!sys.angles.empty()) {
        out << "\nAngles\n\n";
        for (const auto& a : sys.angles) out << a.id << ' ' << a.type << ' ' << a.a << ' ' << a.b << ' ' << a.c << '\n';
    }
    if (!sys.dihedrals.empty()) {
        out << "\nDihedrals\n\n";
        for (const auto& d : sys.dihedrals) out << d.id << ' ' << d.type << ' ' << d.a << ' ' << d.b << ' ' << d.c << ' ' << d.d << '\n';
    }
    if (!out) throw std::runtime_error("Failed while writing output file: " + s.output);
}

} // namespace

int main(int argc, char** argv) {
    try {
        Settings settings = parse_args(argc, argv);
        apply_crosslinker_stoichiometry(settings);
        apply_filler_composition(settings);
        validate(settings);
        report_composition(settings);
        const long long beads = 1LL*settings.n1*settings.m1 + 1LL*settings.n2*settings.m2
                              + 1LL*settings.n3*settings.m3 + 1LL*settings.n4*settings.m4;
        const double box = std::cbrt(beads * settings.bead_mass / (settings.density * kAvogadroScale));
        const System system = build_system(settings, box);
        write_data(settings, system, box);
        std::cerr << "Wrote " << settings.output << " (" << system.atoms.size() << " atoms, "
                  << system.bonds.size() << " bonds, " << system.angles.size() << " angles, "
                  << system.dihedrals.size() << " dihedrals; box " << box << " A)\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
