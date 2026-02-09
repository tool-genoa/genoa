#!/usr/bin/env bash
# =================================================================
# Run GENOA v3 with a given configuration file.
# Captures the output in a log file.
# Optionally sends an email notification upon completion.
#
# Usage:
#   ./runtime/scripts/run_job.sh [config_file] [output_log]
#   ./runtime/scripts/run_job.sh --help
#
# Examples:
#   ./runtime/scripts/run_job.sh
#   ./runtime/scripts/run_job.sh ./test_cases/config/1_build.ini
#   ./runtime/scripts/run_job.sh ./test_cases/config/1_build.ini log.1.txt
#   nohup ./runtime/scripts/run_job.sh ./test_cases/config/1_build.ini log.1.txt &
# =================================================================

# ---- user settings ----
DEFAULT_INPUT="./test_cases/config/1_build.ini"
DEFAULT_OUTPUT="./runtime/logs/log.genoa_output.txt"
GENOA_CMD="python3 -m genoa"    # or GENOA_CMD="genoa"
EMAIL=""                        # set to your_email@domain.com to enable email
# ------------------------

# ---- help ----
if [[ "$1" == "--help" || "$1" == "-h" ]]; then
  echo "Usage: ./runtime/scripts/run_job.sh [config_file] [output_log]"
  echo
  echo "Runs GENOA v3 using the specified configuration file."
  echo
  echo "Arguments:"
  echo "  config_file   Path to the configuration file (default: $DEFAULT_INPUT)"
  echo "  output_log    Log file path (default: $DEFAULT_OUTPUT)"
  echo
  echo "Examples:"
  echo "  ./runtime/scripts/run_job.sh"
  echo "  ./runtime/scripts/run_job.sh ./test_cases/config/1_build.ini log.1.txt"
  echo "  nohup ./runtime/scripts/run_job.sh ./test_cases/config/1_build.ini log.1.txt &"
  exit 0
fi

# ---- args ----
INPUT_FILE="${1:-$DEFAULT_INPUT}"
OUTPUT_FILE="${2:-$DEFAULT_OUTPUT}"

# ---- checks ----
if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Error: configuration file not found: $INPUT_FILE"
  echo "Use './runtime/scripts/run_job.sh --help' for usage information."
  exit 1
fi

mkdir -p "$(dirname "$OUTPUT_FILE")"

# ---- run ----
echo "Running GENOA..."
echo "  Config: $INPUT_FILE"
echo "  Log:    $OUTPUT_FILE"

nohup $GENOA_CMD "$INPUT_FILE" > "$OUTPUT_FILE" 2>&1
status=$?

# ---- status ----
if [[ $status -eq 0 ]]; then
  msg="GENOA completed successfully: $INPUT_FILE"
else
  msg="GENOA encountered an error: $INPUT_FILE (exit code $status)"
fi

# ---- email ----
if [[ -n "$EMAIL" ]] && command -v mail >/dev/null 2>&1; then
  echo "$msg" | mail -s "GENOA job status: $OUTPUT_FILE" "$EMAIL"
fi

echo "$msg"
exit $status
