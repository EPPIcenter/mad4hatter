#!/usr/bin/env python
# coding: utf-8

import pandas as pd


def create_amplicon_info_from_pools(pools, amplicon_info_paths):
    amplicon_info_per_pool = []
    for i in range(len(pools)):
        amplicon_info_for_pool = pd.read_csv(amplicon_info_paths[i])
        amplicon_info_for_pool['pool'] = pools[i]
        amplicon_info_per_pool.append(amplicon_info_for_pool)
    return pd.concat(amplicon_info_per_pool)
