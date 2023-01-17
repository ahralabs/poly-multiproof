#!/bin/bash
RUSTDOCFLAGS="--html-in-header src/docs-header.html" cargo doc --document-private-items
