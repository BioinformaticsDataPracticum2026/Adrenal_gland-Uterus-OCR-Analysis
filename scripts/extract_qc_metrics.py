"""
Extract TSS enrichment, FRiP, and library complexity metrics
from ENCODE ATAC-seq pipeline output directories on Bridges-2.
"""
import json, os

samples = {
    "Human_Adrenal": "/ocean/projects/bio200034p/ikaplow/HumanDNase/AdrenalGlandFemaleGoodReps/atac/eebff493-9ba0-490d-9763-c3339e2de395",
    "Human_Uterus":  "/ocean/projects/bio200034p/ikaplow/HumanDNase/UterusFemale/atac/715c0657-2968-4aad-9a72-e4ca10ff9a34",
    "Mouse_Adrenal": "/ocean/projects/bio200034p/ikaplow/MouseDNase/AdrenalGland2025/atac/f96b2348-7f72-45c7-97f6-65ec93bebb19",
    "Mouse_Uterus":  "/ocean/projects/bio200034p/ikaplow/MouseDNase/Uterus2025/atac/ae91da79-0594-449e-8215-0d10629d55d7",
}

for name, base in samples.items():
    print(f"\n=== {name} ===")
    # TSS enrichment per replicate
    for shard in [0, 1]:
        tss_dir = f"{base}/call-tss_enrich/shard-{shard}/execution"
        if os.path.isdir(tss_dir):
            for f in os.listdir(tss_dir):
                if f.endswith(".tss_enrich.qc"):
                    val = open(f"{tss_dir}/{f}").read().strip()
                    print(f"  TSS Rep{shard+1}: {val}")
    # Pooled FRiP
    frip = f"{base}/call-call_peak_pooled/execution/rep.pooled.pval0.01.300K.bfilt.frip.qc"
    if os.path.exists(frip):
        print(f"  FRiP (pooled): {open(frip).read().strip()}")
    # Library complexity and alignment from qc.json
    qc_path = f"{base}/call-qc_report/execution/qc.json"
    if os.path.exists(qc_path):
        qc = json.load(open(qc_path))
        lc = qc.get("lib_complexity", {})
        for rep in ["rep1", "rep2"]:
            r = lc.get(rep, {})
            print(f"  {rep}: NRF={r.get('NRF','N/A')} PBC1={r.get('PBC1','N/A')} PBC2={r.get('PBC2','N/A')}")
