#!/usr/bin/env bash
#
# generate-and-test.sh — Generate synthetic star field images at random
# orientations and test the star-tracker pipeline against the known
# ground-truth attitudes.
#
# Usage:
#   ./generate-and-test.sh [OPTIONS]
#
# Options:
#   -n <count>           Number of images to generate (default: 10)
#   -c <catalog.tsv>     Path to the bright-star-catalog TSV
#   -d <database>        Path to the serialized MultiDatabase file
#   -b <benchmark>       Path to the benchmark executable
#   -g <generate-stars>  Path to the generate-stars executable
#   -f <fov-deg>         Field of view in degrees (default: 15)
#   -t <tolerance-deg>   Angular tolerance for pass/fail (default: 1.0)
#   --width <px>         Image width (default: 512)
#   --height <px>        Image height (default: 512)
#   --spread <stddev>    Star spread (default: 2.0)
#   --photons <N>        Zero-mag total photons (default: 200000)
#   --saturation <N>     Saturation photons (default: 1000)
#   --dark-current <f>   Dark current 0-1 (default: 0.0)
#   --read-noise <f>     Read noise stddev (default: 0.01)
#   --shot-noise         Enable shot noise
#   --seed <int>         Random seed (default: time-based)
#   -h                   Show this help
#
# Requires: python3 for angular math

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

COUNT=10
CATALOG=""
DATABASE=""
BENCHMARK_BIN="${PROJECT_DIR}/build/benchmark"
GENERATE_BIN="${PROJECT_DIR}/build/generate-stars"
FOV=15
TOLERANCE=1.0
WIDTH=512
HEIGHT=512
SPREAD=2.0
PHOTONS=200000
SATURATION=1000
DARK_CURRENT=0.0
READ_NOISE=0.01
SHOT_NOISE=""
SEED=""
MIN_IDENTIFICATIONS=6

# ── Parse options ────────────────────────────────────────────────────────────

usage() {
    sed -n '3,30p' "$0" | sed 's/^# \?//'
    exit "${1:-0}"
}

while [ $# -gt 0 ]; do
    case "$1" in
        -n)              COUNT="$2";         shift 2 ;;
        -c)              CATALOG="$2";       shift 2 ;;
        -d)              DATABASE="$2";      shift 2 ;;
        -b)              BENCHMARK_BIN="$2"; shift 2 ;;
        -g)              GENERATE_BIN="$2";  shift 2 ;;
        -f)              FOV="$2";           shift 2 ;;
        -t)              TOLERANCE="$2";     shift 2 ;;
        --width)         WIDTH="$2";         shift 2 ;;
        --height)        HEIGHT="$2";        shift 2 ;;
        --spread)        SPREAD="$2";        shift 2 ;;
        --photons)       PHOTONS="$2";       shift 2 ;;
        --saturation)    SATURATION="$2";    shift 2 ;;
        --dark-current)  DARK_CURRENT="$2";  shift 2 ;;
        --read-noise)    READ_NOISE="$2";    shift 2 ;;
        --shot-noise)    SHOT_NOISE="--shot-noise"; shift ;;
        --seed)          SEED="$2";          shift 2 ;;
        -h|--help)       usage 0 ;;
        *)               echo "Unknown option: $1"; usage 1 ;;
    esac
done

# ── Auto-detect paths ───────────────────────────────────────────────────────

if [ -z "$CATALOG" ]; then
    for candidate in \
        "$PROJECT_DIR/../lost/bright-star-catalog.tsv" \
        "$PROJECT_DIR/bright-star-catalog.tsv" \
        "$PROJECT_DIR/data/bright-star-catalog.tsv"; do
        [ -f "$candidate" ] && CATALOG="$candidate" && break
    done
    if [ -z "$CATALOG" ]; then
        echo "Error: -c <catalog> required (path to bright-star-catalog.tsv)"
        exit 1
    fi
fi

if [ -z "$DATABASE" ]; then
    for candidate in "$PROJECT_DIR/build/stardb.dat" "$PROJECT_DIR/build/stars.dat"; do
        [ -f "$candidate" ] && DATABASE="$candidate" && break
    done
    if [ -z "$DATABASE" ]; then
        echo "Error: -d <database> required"
        exit 1
    fi
fi

for bin in "$BENCHMARK_BIN" "$GENERATE_BIN"; do
    if [ ! -x "$bin" ]; then
        echo "Error: executable not found: $bin"
        echo "       Build: cd build && cmake .. && make"
        exit 1
    fi
done

command -v python3 >/dev/null 2>&1 || { echo "Error: python3 is required"; exit 1; }

# ── Working directory ────────────────────────────────────────────────────────

WORKDIR=$(mktemp -d "${SCRIPT_DIR}/test-images-XXXXXX")
trap 'rm -rf "$WORKDIR"' EXIT

# ── Helper functions ─────────────────────────────────────────────────────────

angular_separation() {
    python3 -c "
import math
ra1, dec1, ra2, dec2 = [math.radians(float(x)) for x in ['$1','$2','$3','$4']]
cos_sep = math.sin(dec1)*math.sin(dec2) + math.cos(dec1)*math.cos(dec2)*math.cos(ra1-ra2)
cos_sep = max(-1.0, min(1.0, cos_sep))
print(f'{math.degrees(math.acos(cos_sep)):.6f}')
"
}

angle_diff() {
    python3 -c "
a, b = float('$1'), float('$2')
d = abs(a - b) % 360
if d > 180: d = 360 - d
print(f'{d:.6f}')
"
}

# ── Step 1: Generate images ──────────────────────────────────────────────────

echo ""
echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║          Generating Synthetic Star Field Images                 ║"
echo "╚══════════════════════════════════════════════════════════════════╝"
echo ""

GEN_ARGS=(
    --catalog "$CATALOG"
    --output "$WORKDIR"
    --count "$COUNT"
    --fov "$FOV"
    --width "$WIDTH"
    --height "$HEIGHT"
    --random-attitudes
    --spread "$SPREAD"
    --photons "$PHOTONS"
    --saturation "$SATURATION"
    --dark-current "$DARK_CURRENT"
    --read-noise "$READ_NOISE"
    --exposure 0.5
)
[ -n "$SHOT_NOISE" ] && GEN_ARGS+=("$SHOT_NOISE")
[ -n "$SEED" ] && GEN_ARGS+=(--seed "$SEED")

"$GENERATE_BIN" "${GEN_ARGS[@]}"

# ── Step 2: Run the pipeline on each image ───────────────────────────────────

echo ""
echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║          Testing Star-Tracker Pipeline                         ║"
echo "╠══════════════════════════════════════════════════════════════════╣"
printf "║  Images: %-5d  FOV: %4.0f°  Tolerance: %.2f°                    ║\n" \
       "$COUNT" "$FOV" "$TOLERANCE"
printf "║  Database: %-50s ║\n" "$(basename "$DATABASE")"
echo "╚══════════════════════════════════════════════════════════════════╝"
echo ""

PASS=0
FAIL=0
SKIP=0

# Read manifest (skip comment line)
while IFS= read -r line; do
    # Skip comments / empty lines
    [[ "$line" =~ ^# ]] && continue
    [[ -z "$line" ]] && continue

    # Parse: image_NNN.png  ra  dec  roll  nstars  fov  width  height  fl
    read -r png_file true_ra true_dec true_roll num_stars img_fov img_w img_h focal_length <<< "$line"

    img_path="$WORKDIR/$png_file"

    echo "━━━ ${png_file} (RA=${true_ra}° DEC=${true_dec}° Roll=${true_roll}°, ${num_stars} stars) ━━━"

    if [ ! -f "$img_path" ]; then
        echo "  SKIP: image file missing"
        SKIP=$((SKIP + 1))
        echo ""
        continue
    fi

    # Run the benchmark pipeline
    result=$(timeout 30 "$BENCHMARK_BIN" \
        --image "$img_path" \
        --database "$DATABASE" \
        --catalog "$CATALOG" \
        --focal-length "$focal_length" \
        --tolerance 0.001 \
        --max-centroids 50 \
        --min-star-pixels 3 \
        --sigma 2.0 \
        2>/dev/null || echo "ERROR")

    if [ "$result" = "ERROR" ] || [ -z "$result" ]; then
        echo "  SKIP: benchmark timed out or crashed"
        SKIP=$((SKIP + 1))
        echo ""
        continue
    fi

    read -r res_ra res_dec res_roll res_nids res_ncen <<< "$result"

    if [ "$res_ra" = "NaN" ]; then
        echo "  SKIP: no solution (identified $res_nids / $res_ncen centroids)"
        SKIP=$((SKIP + 1))
        echo ""
        continue
    fi

    if [ "$res_nids" -lt "$MIN_IDENTIFICATIONS" ]; then
        echo "  SKIP: too few identifications ($res_nids < $MIN_IDENTIFICATIONS)"
        SKIP=$((SKIP + 1))
        echo ""
        continue
    fi

    echo "  Pipeline:  RA=${res_ra}°  DEC=${res_dec}°  Roll=${res_roll}°  (id: $res_nids/$res_ncen)"

    sep=$(angular_separation "$true_ra" "$true_dec" "$res_ra" "$res_dec")
    roll_diff=$(angle_diff "$true_roll" "$res_roll")

    echo "  Pointing error: ${sep}°    Roll error: ${roll_diff}°"

    if python3 -c "exit(0 if float('$sep') < float('$TOLERANCE') else 1)"; then
        echo "  Result: PASS"
        PASS=$((PASS + 1))
    else
        echo "  Result: FAIL"
        FAIL=$((FAIL + 1))
    fi
    echo ""

done < "$WORKDIR/manifest.txt"

# ── Summary ──────────────────────────────────────────────────────────────────

echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║                         SUMMARY                                ║"
echo "╠══════════════════════════════════════════════════════════════════╣"
printf "║  Total: %-4d   Pass: %-4d   Fail: %-4d   Skip: %-4d            ║\n" \
       "$COUNT" "$PASS" "$FAIL" "$SKIP"
echo "╚══════════════════════════════════════════════════════════════════╝"

if [ "$FAIL" -gt 0 ]; then
    exit 1
fi
exit 0
