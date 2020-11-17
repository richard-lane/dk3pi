#include "cuts.h"

int main()
{
    semileptonicCuts("wg_rs_sl.root", "DecayTree", "cut_wg_rs_sl.root", "D0_PT", "D0_M");

    promptCuts("wg_rs_prompt.root",
               "DecayTree",
               "cut_wg_rs_prompt.root",
               "Dst_pi_PT",
               "Dst_pi_TRACK_GhostProb",
               "D0_PT",
               "D0_M");
}
