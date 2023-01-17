#!/bin/bash
cargo criterion --plotting-backend disabled -- --quick --quiet &> bench_out.txt
