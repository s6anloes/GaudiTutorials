// k4FWCore
#include "k4FWCore/Consumer.h"

// edm4hep
#include "edm4hep/SimCalorimeterHitCollection.h"

// STL
#include <string>
#include <vector>
#include <cmath>

struct MoliereRadiusSolution final
    : k4FWCore::Consumer<void(const edm4hep::SimCalorimeterHitCollection&)> {

public:
    // Constructor
    MoliereRadiusSolution(const std::string& name, ISvcLocator* svcLoc)
        : Consumer(name, 
                   svcLoc, 
                   KeyValues("InputCaloHitCollection", {"simplecaloRO"})) {}

    
    void operator()(const edm4hep::SimCalorimeterHitCollection& input) const override {

        // Calcualte the Moliere radius from the simhits

        double barycentre_x = 0.0;
        double barycentre_y = 0.0;
        double totalEnergy = 0.0;

        for (const auto& simhit: input) {
            // Calculate the shower barycentre in x and y
            double x = simhit.getPosition().x;
            double y = simhit.getPosition().y;

            double energy = simhit.getEnergy();

            barycentre_x += x * energy;
            barycentre_y += y * energy;
            totalEnergy += energy;
        }

        if (totalEnergy == 0.0) {
            info() << "Total energy is zero; cannot compute Moliere radius." << endmsg;
            return;
        }

        barycentre_x /= totalEnergy;
        barycentre_y /= totalEnergy;

        std::vector<HitData> radialHits;

        for (const auto& hit : input) {
            const auto& pos = hit.getPosition();
            double dx = pos.x - barycentre_x;
            double dy = pos.y - barycentre_y;
            double r = std::sqrt(dx * dx + dy * dy);
            double e = hit.getEnergy();
            radialHits.push_back({r, e});
        }

        // Step 3: Sort hits by radius
        std::sort(radialHits.begin(), radialHits.end(), [](const HitData& a, const HitData& b) {
            return a.r < b.r;
        });

        // Step 4: Accumulate energy and find radius at 90%
        double cumulativeEnergy = 0.0;
        double moliereRadius = 0.0;
        for (const auto& hit : radialHits) {
            cumulativeEnergy += hit.energy;
            if (cumulativeEnergy >= 0.9 * totalEnergy) {
                moliereRadius = hit.r;
                break;
            }
        }

        info() << "Moliere radius = " << moliereRadius << " mm" << endmsg;

        
    }

private:
    struct HitData {
        double r;      // Radius from barycentre
        double energy; // Energy of the hit
    };
};

DECLARE_COMPONENT(MoliereRadiusSolution)