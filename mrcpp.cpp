#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "MRCPP/Plotter"

#include "MRCPP/MWFunctions"
#include "MRCPP/Gaussians"
#include "MRCPP/MWOperators"

#include <fstream>

// This is your personal sandbox to test out ideas.
// See the examples/ directory for inspiration.
// Please do not commit your pilot code to git.

using namespace mrcpp;
const auto min_scale = 0;
const auto max_depth = 25;

const auto order = 7;
const auto prec = 1.0e-5;
const auto D = 3; // Dimensions
int main(int argc, char **argv) {
    Timer timer;

    // Initialize printing
    int printlevel = 0;
    Printer::init(printlevel);
    Printer::printEnvironment();
    Printer::printHeader(0, "MRCPP pilot code");

    // Plots
    int nPts = 1000;
    double a[3] = {0.0, 0.0, 0.0};
    double b[3] = {20.0, 0.0, 0.0};

    Plotter<3> plot;
    plot.setNPoints(nPts);
    plot.setRange(a, b);

{
    auto basis = InterpolatingBasis(order);
    auto sf_x = 2.0*pi;
    auto sf = std::array<double, D>{sf_x, sf_x, sf_x};
    auto world = BoundingBox<D>(sf, true);

    // auto world = BoundingBox<D>(sf, true);
    auto MRA = MultiResolutionAnalysis<D>(world, basis, max_depth);


    auto g_func = [] (const mrcpp::Coord<D> &r) {
        // auto p = pi*1.0;
        // return cos(p*r[0])*cos(p*r[1])*cos(p*r[2]);
        return cos(r[0])*cos(r[1])*cos(r[2]);
    };
    auto dg_func = [] (const mrcpp::Coord<D> &r) {
        // auto p = pi*1.0;
        // return cos(p*r[0])*cos(p*r[1])*cos(p*r[2]);
        return sin(r[0])*cos(r[1])*cos(r[2]);
    };

    FunctionTree<D> g_tree(MRA);
    FunctionTree<D> dg_tree(MRA);
    FunctionTree<D> dg_exact_tree(MRA);
    ABGVOperator<D> Der(MRA, 0.0, 0.0);

    project<D>(prec, g_tree, g_func);
    project<D>(prec, dg_exact_tree, dg_func);

   apply(dg_tree, Der, g_tree, 0);


    // println(0, f_tree.evalf({0.0, 0.0, 0.0}));
    println(0, g_tree.evalf({0.0, 0.0, 0.0}));
    plot.linePlot(dg_tree, "dg_tree_sf");
    MRA.print();
}
// {
//     auto basis = InterpolatingBasis(order);
//     auto sf_x = 1.0;
//     auto sf = std::array<double, D>{sf_x, sf_x, sf_x};
//     auto world = BoundingBox<D>(-4, {0, 0, 0}, {2, 2, 2});
//
//     // auto world = BoundingBox<D>(sf, true);
//     auto MRA = MultiResolutionAnalysis<D>(world, basis, max_depth);
//
//     auto f_func = [] (const mrcpp::Coord<D> &r) {
//         auto p = 1.0;
//         return cos(p*r[0])*cos(p*r[1])*cos(p*r[2]);
//     };
//
//     FunctionTree<D> f_tree(MRA);
//     FunctionTree<D> g_tree(MRA);
//     PoissonOperator P(MRA, prec);
//     project<D>(prec, f_tree, f_func);
//
//     apply(prec, g_tree, P, f_tree);
//     plot.linePlot(g_tree, "g_tree");
//
//     MRA.print();
// }
    timer.stop();
    Printer::printFooter(0, timer, 2);
    return 0;
}
