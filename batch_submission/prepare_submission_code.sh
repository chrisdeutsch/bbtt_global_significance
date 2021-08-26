#!/usr/bin/env bash
tar --exclude-vcs --exclude='*.pcm' --exclude='*.d' --exclude='*.so' \
    -czvf /cephfs/user/s6crdeut/bbtt_global_significance.tar.gz \
    -C ../../ bbtt_global_significance/{scripts,macros}
