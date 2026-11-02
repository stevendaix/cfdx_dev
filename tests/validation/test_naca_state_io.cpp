#include "naca_qualification_state_io.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
using cfdx::core::Vec3;
using cfdx::validation::NacaStateRow;
using cfdx::validation::naca_read_state_csv;
using cfdx::validation::naca_write_state_csv;
static void require(bool x) { if (!x) throw std::runtime_error("state I/O contract failed"); }
int main() {
    const std::string p = "/tmp/cfdx_naca_state_io_contract.csv";
    const std::vector<NacaStateRow> states = {
        {{0.0,0.0,0.5},1.0,0.0,-0.03},
        {{1.0,0.1,0.5},0.3,0.02,0.4}};
    naca_write_state_csv(p, states);
    auto load = naca_read_state_csv(p, {states[0].center, states[1].center});
    require(load.size() == states.size());
    require(std::abs(load[1].pressure - 0.4) < 1e-15);
    // A Windows CSV writer may use CRLF: import must not fail spuriously.
    const std::string p_crlf = p + ".crlf";
    {
        std::ifstream source(p, std::ios::binary);
        std::ofstream target(p_crlf, std::ios::binary);
        char ch;
        while (source.get(ch)) {
            if (ch == '\n') target.put('\r');
            target.put(ch);
        }
    }
    load = naca_read_state_csv(p_crlf, {states[0].center, states[1].center});
    require(load.size() == 2 && std::abs(load[0].ux - 1.0) < 1e-15);
    std::remove(p_crlf.c_str());
    const std::string atom = p + ".atomic";
    cfdx::validation::naca_write_state_csv_atomic(atom, states);
    require(naca_read_state_csv(atom, {states[0].center, states[1].center}).size() == 2);
    { std::ifstream incomplete(atom + ".incomplete"); require(!incomplete); }
    std::remove(atom.c_str());
    bool fail = false;
    try { (void)naca_read_state_csv(p, {{0,0,0.5},{1,0.2,0.5}}); }
    catch (const std::invalid_argument&) { fail = true; }
    require(fail);
    fail = false;
    try { (void)naca_read_state_csv(p, {states[0].center}); }
    catch (const std::invalid_argument&) { fail = true; }
    require(fail);
    std::remove(p.c_str());
    std::cout << "NACA_STATE_IO: PASS\n";
}
