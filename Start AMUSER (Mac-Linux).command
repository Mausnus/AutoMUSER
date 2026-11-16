#!/bin/bash
# ====================================================================
#  Start the AMUSER GUI on macOS or Linux.
#  macOS: double-click this file (you may need to right-click > Open
#         the first time, to get past Gatekeeper).
#  Linux: run  ./"Start AMUSER (Mac-Linux).command"  in a terminal,
#         or double-click if your file manager allows it.
#  A browser tab opens automatically. Close this window to stop.
# ====================================================================
cd "$(dirname "$0")" || exit 1
if ! command -v perl >/dev/null 2>&1; then
  echo "Perl was not found. macOS ships with Perl; on Linux install it with"
  echo "your package manager (e.g. sudo apt install perl), then try again."
  read -r -p "Press Enter to close."
  exit 1
fi
perl amuser_gui.pl
echo
echo "The AMUSER server has stopped."
read -r -p "Press Enter to close."
