#include "DQCD/Modules/interface/EffGen.h"

EffGen::EffGen (std::string filename, std::string tag):
    corr(filename, tag)
{}

// Destructor
EffGen::~EffGen() {}

float EffGen::get_weight(
    int displaced_pdgid, float input_ctau, float output_ctau,
    iRVec GenPart_pdgId, iRVec GenPart_genPartIdxMother, fRVec GenPart_mass,
    fRVec GenPart_pt, fRVec GenPart_vertex_x, fRVec GenPart_vertex_y, fRVec GenPart_vertex_z,
    fRVec GenPart_eta, fRVec GenPart_phi
) {
    float num_weight = 1.;
    float den_weight = 1.;
    std::vector<int> indexes;
    for (size_t i = 0; i < GenPart_pdgId.size(); i++) {
        if (abs(GenPart_pdgId[i]) != 13) {
            continue;
        }
        auto mother_index = GenPart_genPartIdxMother[i];
        if (mother_index == -1)
            continue;
        if (std::find(indexes.begin(), indexes.end(), mother_index) != std::end(indexes)) {
            // Avoid considering the same LLP more than once
            continue;
        }
        if (GenPart_pt[mother_index] < 5)
            continue;
        // Is the mother particle the LLP we are looking for?
        if (abs(GenPart_pdgId[mother_index]) == displaced_pdgid) {
            indexes.push_back(mother_index);

            auto lxy = std::sqrt(
                std::pow(GenPart_vertex_x[i] - GenPart_vertex_x[mother_index], 2) +
                std::pow(GenPart_vertex_y[i] - GenPart_vertex_y[mother_index], 2)
            );
			auto lxyz = std::sqrt(
                std::pow(GenPart_vertex_x[i] - GenPart_vertex_x[mother_index], 2) +
                std::pow(GenPart_vertex_y[i] - GenPart_vertex_y[mother_index], 2) +
                std::pow(GenPart_vertex_z[i] - GenPart_vertex_z[mother_index], 2) 
            );

            //auto px = GenPart_pt[mother_index] * cos(GenPart_phi[mother_index]);
            //auto py = GenPart_pt[mother_index] * sin(GenPart_phi[mother_index]);
            //auto pz = GenPart_pt[mother_index] * sinh(GenPart_eta[mother_index]);

            auto p =  GenPart_pt[mother_index] * cosh(GenPart_eta[mother_index]);

            // auto lxy_pt = (
                // (GenPart_vertex_x[i] - GenPart_vertex_x[mother_index]) * px +
                // (GenPart_vertex_y[i] - GenPart_vertex_y[mother_index]) * py
            // );
            // auto ct = lxy_pt * GenPart_mass[mother_index] / std::pow(GenPart_pt[mother_index], 2);
            auto ct = 10 * GenPart_mass[mother_index] * lxyz / p;
            // auto ct = lxyz * GenPart_mass[mother_index] / std::sqrt(px * px + py * py + pz * pz) ;
            // auto ct = lxy / (GenPart_mass[mother_index] * GenPart_pt[mother_index]);
			// std::cout << ct << " " << lxy << " ";
            num_weight *= 1 - ((corr.eval({lxy, GenPart_pt[mother_index], "sf"})) * (
                (1. / output_ctau) * std::exp(-ct / output_ctau)
            ));
            den_weight *= 1 - ((corr.eval({lxy, GenPart_pt[mother_index], "sf"})) * (
                (1. / input_ctau) * std::exp(-ct / input_ctau)
            ));
			// std::cout << num_weight << " " << den_weight << std::endl;
        }
    }
    if (den_weight == 1.)
        return 0.;
    else
        return (1 - num_weight) / (1 - den_weight);
}


float EffGen::get_weight_noeff(
    int displaced_pdgid, float input_ctau, float output_ctau,
    iRVec GenPart_pdgId, iRVec GenPart_genPartIdxMother, fRVec GenPart_mass,
    fRVec GenPart_pt, fRVec GenPart_vertex_x, fRVec GenPart_vertex_y, fRVec GenPart_vertex_z,
    fRVec GenPart_eta, fRVec GenPart_phi
) {
    float total_weight = 0.;
    size_t nA = 0;
    std::vector<int> indexes;
    for (size_t i = 0; i < GenPart_pdgId.size(); i++) {
        if (abs(GenPart_pdgId[i]) != 13) {
            continue;
        }
        auto mother_index = GenPart_genPartIdxMother[i];
        if (mother_index == -1)
            continue;
        if (std::find(indexes.begin(), indexes.end(), mother_index) != std::end(indexes)) {
            // Avoid considering the same LLP more than once
            continue;
        }
        if (GenPart_pt[mother_index] < 15)
            continue;
        // Is the mother particle the LLP we are looking for?
        if (abs(GenPart_pdgId[mother_index]) == displaced_pdgid) {
            indexes.push_back(mother_index);

            auto lxy = std::sqrt(
                std::pow(GenPart_vertex_x[i] - GenPart_vertex_x[mother_index], 2) +
                std::pow(GenPart_vertex_y[i] - GenPart_vertex_y[mother_index], 2)
            );
			auto lxyz = std::sqrt(
                std::pow(GenPart_vertex_x[i] - GenPart_vertex_x[mother_index], 2) +
                std::pow(GenPart_vertex_y[i] - GenPart_vertex_y[mother_index], 2) +
                std::pow(GenPart_vertex_z[i] - GenPart_vertex_z[mother_index], 2) 
            );

            //auto px = GenPart_pt[mother_index] * cos(GenPart_phi[mother_index]);
            //auto py = GenPart_pt[mother_index] * sin(GenPart_phi[mother_index]);
            //auto pz = GenPart_pt[mother_index] * sinh(GenPart_eta[mother_index]);

            auto p =  GenPart_pt[mother_index] * cosh(GenPart_eta[mother_index]);

            // auto lxy_pt = (
                // (GenPart_vertex_x[i] - GenPart_vertex_x[mother_index]) * px +
                // (GenPart_vertex_y[i] - GenPart_vertex_y[mother_index]) * py
            // );
            // auto ct = lxy_pt * GenPart_mass[mother_index] / std::pow(GenPart_pt[mother_index], 2);
            auto ct = 10 * GenPart_mass[mother_index] * lxyz / p;
            // auto ct = lxyz * GenPart_mass[mother_index] / std::sqrt(px * px + py * py + pz * pz) ;
            // auto ct = lxy / (GenPart_mass[mother_index] * GenPart_pt[mother_index]);
			// std::cout << ct << " " << lxy << " ";

            nA++;
            total_weight = ((1. / output_ctau) * std::exp(-ct / output_ctau))
                / ((1. / input_ctau) * std::exp(-ct / input_ctau));
        }
    }
    if (nA == 0)
        return total_weight;  // 0
    else
        return total_weight / float(nA);
}

