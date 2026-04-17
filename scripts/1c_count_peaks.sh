#!/bin/bash
# Count IDR peaks for all four ATAC-seq samples
HUMAN_ADRENAL="/ocean/projects/bio200034p/ikaplow/HumanDNase/AdrenalGlandFemaleGoodReps/atac/eebff493-9ba0-490d-9763-c3339e2de395/call-reproducibility_idr/execution/glob-aae3f88ec555ee05b83e37921e9eb318/idr.conservative_peak.narrowPeak.gz"
HUMAN_UTERUS="/ocean/projects/bio200034p/ikaplow/HumanDNase/UterusFemale/atac/715c0657-2968-4aad-9a72-e4ca10ff9a34/call-reproducibility_idr/execution/glob-aae3f88ec555ee05b83e37921e9eb318/idr.conservative_peak.narrowPeak.gz"
MOUSE_ADRENAL="/ocean/projects/bio200034p/ikaplow/MouseDNase/AdrenalGland2025/atac/f96b2348-7f72-45c7-97f6-65ec93bebb19/call-reproducibility_idr/execution/glob-aae3f88ec555ee05b83e37921e9eb318/idr.conservative_peak.narrowPeak.gz"
MOUSE_UTERUS="/ocean/projects/bio200034p/ikaplow/MouseDNase/Uterus2025/atac/ae91da79-0594-449e-8215-0d10629d55d7/call-reproducibility_idr/execution/glob-aae3f88ec555ee05b83e37921e9eb318/idr.conservative_peak.narrowPeak.gz"

for f in "$HUMAN_ADRENAL" "$HUMAN_UTERUS" "$MOUSE_ADRENAL" "$MOUSE_UTERUS"; do
    echo "$f: $(zcat $f | wc -l) peaks"
done
