#!/usr/bin/env python3
"""WfAverage — CD waveform event-averaging runner."""

import argparse
import Sniper


def get_parser():
    parser = argparse.ArgumentParser(description="CD waveform event-average")
    parser.add_argument("--evtmax", type=int, default=-1,
                        help="events to process (-1 = all)")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--input", nargs="+",
                       help="input rtraw file(s)")
    group.add_argument("--input-list",
                       help="text file listing input rtraw paths, one per line")
    parser.add_argument("--user-output", default="wf_avg.root",
                        help="output ROOT file")
    return parser


def run(args):
    if args.input_list:
        with open(args.input_list) as f:
            input_files = [line.strip() for line in f if line.strip()]
    else:
        input_files = args.input

    task = Sniper.Task("WfAverageTask")
    task.setEvtMax(args.evtmax)
    task.setLogLevel(2)

    import BufferMemMgr
    bufMgr = task.createSvc("BufferMemMgr")
    bufMgr.property("TimeWindow").set([0, 0])

    import RootIOSvc
    inputsvc = task.createSvc("RootInputSvc/InputSvc")
    inputsvc.property("InputFile").set(input_files)

    import RootWriter
    rootwriter = task.createSvc("RootWriter")
    rootwriter.property("Output").set({"USER_OUTPUT": args.user_output})

    import Geometry
    geom = task.createSvc("RecGeomSvc")
    geom.property("GeomFile").set("default")
    geom.property("FastInit").set(True)

    Sniper.loadDll("build/lib/libWfAverage.so")
    task.createAlg("WfAverage")

    task.show()
    task.run()
    del task


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    run(args)
