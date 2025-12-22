#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys
from diya import (diya, load_power_matrix_from_csv, prepare_protocol,
                  load_device_config, generate_color_help_text, NROWS, NCOLS,
                  get_reference_powers, load_device_calibration,
                  upload_protocol)

##############################################################################
########################### USER CONFIGURATION ###############################
##############################################################################
# These are the default values that will be used if not specified via command line

# Optical power matrix, in mW/cm2 (if calibrated) 3D array: [channel, row, column]
# Note: Colors 0-3 use channels 0-3 with LEDG1, Colors 4-7 use channels 0-3 with LEDG2
DEFAULT_OPTICAL_POWER = np.array([
	# Channel 0 (Color 0 or 4). Blue on v0.4c
	[[0 / (2**(7-row)) for col in range(NCOLS)] for row in range(NROWS)],
	# Channel 1 (Color 1 or 5). Green or Violet on v0.4c
	np.array([[x*row for x in [0,0,0,0,0,0,0.3,1,1.5,2,2.8,0]] for row in [0,1,1,1,1,1,1,1]]),
	# Channel 2 (Color 2 or 6). Yellow-Green or White on v0.4c
	[[0 / (2**row) for col in range(NCOLS)] for row in range(NROWS)],
	# Channel 3 (Color 3 or 7). Red or Infrared on v0.4c
	np.array([[x*row for x in [0,2.8,2,1.5,1,0.3,0,0,0,0,0,0]] for row in [0,1,1,1,1,1,1,1]])
])

import numpy as np

DEFAULT_OPTICAL_POWER = np.array([

	# Channel 0 (Color 0 or 4). Blue on v0.4c
	[[0 / (2**(7-row)) for col in range(NCOLS)] for row in range(NROWS)],

    # Channel 1 – Green light (checkerboard starting green at B2)
    np.array([
        # Columns 1–12 (A–L), Rows A–H
        [2.8,	0.56,   1.4,    2.8,	0.56,   1.4,    2.8,	0.56,   1.4,   	2.8,	0.56,   1.4,],  # Row A
        [0.0,  	0.28,  	0.28,  	0.0,  	0.28,  	0.28,	0.0,  	0.28,  	0.28,	0.0,  	0.28,  	0.28,],  # Row B
        [2.8,  	1.4,  	0.56,  	2.8,  	1.4,  	0.56,  	2.8,  	1.4,  	0.56,  	2.8,  	1.4,  	0.56,],  # Row C
        [0.028,	0.028,  2.8,    0.028,	0.028,  2.8,    0.028,	0.028,  2.8,   	0.028,	0.028,  2.8,],  # Row D
        [1.4,  	2.8,  	0.0,  	1.4,  	2.8,  	0.0, 	1.4,  	2.8,  	0.0,    1.4,  	2.8,  	0.0,],  # Row E
        [0.28,	0.0,    2.8,    0.28,	0.0,    2.8,  	0.28,	0.0,    2.8,  	0.28,	0.0,    2.8,],  # Row F
        [0.56,	2.8,    0.028,  0.56,	2.8,    0.028, 	0.56,	2.8,    0.028,  0.56,	2.8,    0.028,],  # Row Gh
        [0.0,	0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # Row H
    ]),
	# Channel 2 (Color 2 or 6). Yellow-Green or White on v0.4c
	[[0 / (2**row) for col in range(NCOLS)] for row in range(NROWS)],    

    # Channel 3 – Red light (opposite checkerboard cells)
    np.array([
        [0.0,	2.8,	2.8,  	0.0,  	2.8,  	2.8,  	0.0,  	2.8,  	2.8,  	0.0,  	2.8,    2.8],  # Row A
        [2.8,   2.8,    2.8,    2.8,    2.8,    2.8,    2.8,    2.8,    2.8,   	2.8,    2.8,    2.8],  # Row B
        [2.8,   2.8,    2.8,    2.8,    2.8,    2.8,    2.8,    2.8,    2.8,   	2.8,    2.8,    2.8],  # Row C
        [2.8,  	2.8,  	0.0,  	2.8,  	2.8,  	0.0,  	2.8,  	2.8,  	0.0,  	2.8,  	2.8,    0.0],  # Row D
        [2.8,   2.8,    2.8,    2.8,   	2.8,    2.8,  	2.8,  	2.8,  	2.8,  	2.8,  	2.8,    2.8],  # Row E
        [2.8,  	2.8,  	2.8,  	2.8,  	2.8,  	2.8, 	2.8,    2.8,    2.8,    2.8,    2.8,    2.8],  # Row F
        [2.8,   0.0,    2.8,    2.8,    0.0,    2.8,    2.8,   	0.0,   	2.8,    2.8,    0.0,    2.8],  # Row G
        [0.0,   0.0,    0.0,   	0.0, 	0.0,  	0.0, 	0.0,  	0.0,  	0.0,  	0.0,  	0.0,   	0.0],  # Row H
    ])
])


print(DEFAULT_OPTICAL_POWER)

""" DEFAULT_OPTICAL_POWER = np.array([
	# Channel 0 (Color 0 or 4). Blue on v0.4c
	[[0 / (2**(7-row)) for col in range(NCOLS)] for row in range(NROWS)],
	# Channel 1 (Color 1 or 5). Green or Violet on v0.4c
	[[2 / (2**col) for col in range(NCOLS)] for row in range(NROWS)],
	# Channel 2 (Color 2 or 6). Yellow-Green or White on v0.4c
	[[0 / (2**row) for col in range(NCOLS)] for row in range(NROWS)],
	# Channel 3 (Color 3 or 7). Red or Infrared on v0.4c
	[[2 / (2**(11-col)) for col in range(NCOLS)] for row in range(NROWS)]
]) """

# Timing parameters (in seconds)
# Must be a multiple of 1/time_multiplier (1/8 = 0.125)
DEFAULT_ON_TIME = 5  # On time in seconds
DEFAULT_PERIOD = 5  # Period in seconds (On time + Off time)

# Temperature setpoint (in °C)
DEFAULT_TEMPERATURE = 37.0

##############################################################################
############################ ADVANCED OPTIONS ################################
##############################################################################
# Default start time matrix (in seconds)
# If not specified, all wells start at the same time (0)
DEFAULT_START_TIME = np.zeros((NROWS, NCOLS))

# Default individual on time matrix (in seconds)
# If not specified, all wells use the same on time from DEFAULT_ON_TIME
DEFAULT_ON_TIME_MATRIX = None  # Will be initialized based on DEFAULT_ON_TIME

# LED group (1 or 2, or both)
DEFAULT_GROUP = 1  # Default to LEDG1 (colors 0-3)

# LED color (0-7). Check hardware led sequence on the top board, left side.
DEFAULT_COLOR = 0  # Default color (legacy mode only)

##############################################################################
######################## END USER CONFIGURATION #############################
##############################################################################

# Add argument parser
parser = argparse.ArgumentParser(description='Upload data to DIYA device with multi-LED-group support')
parser.add_argument('device_id', type=str, help='Device ID (e.g., diya00)')

# Legacy single-colour arguments (for backward compatibility)
parser.add_argument('--power-matrix-csv', type=str, help='[Legacy] Path to CSV file for single color (use with --color)')
parser.add_argument('--start-time-csv', type=str, help='[Legacy] Path to CSV file for start time matrix (use with --color)')
parser.add_argument('--on-time-csv', type=str, help='[Legacy] Path to CSV file for on time matrix (use with --color)')
parser.add_argument('--color', type=int, default=DEFAULT_COLOR, help=f'[Legacy] Single color to use with --power-matrix-csv. {generate_color_help_text()}')

# LED group selection
parser.add_argument('--led-group', type=int, default=DEFAULT_GROUP, choices=[1, 2], help='LED group to enable: 1 (LEDG1, colors 0-3) or 2 (LEDG2, colors 4-7). Default: 1')

# Shared timing and device parameters
parser.add_argument('--on-time', type=float, default=DEFAULT_ON_TIME, help=f'Global on time in seconds (timing resolution auto-calculated or set by --lpwm-value)')
parser.add_argument('--period', type=float, default=DEFAULT_PERIOD, help=f'Period in seconds (timing resolution auto-calculated or set by --lpwm-value)')
parser.add_argument('--lpwm-value', type=int, help='PWM resolution (0-15). If not specified, optimal value will be calculated automatically. Higher values give finer resolution: 0=1s, 3=0.125s, 10=1ms, 15=0.03ms')
parser.add_argument('--temperature', type=float, default=DEFAULT_TEMPERATURE, help=f'Temperature setpoint in Celsius')

# Control flags
parser.add_argument('--no-verify', action='store_true', help='Disable read-back verification (not recommended)')
parser.add_argument('--debug', action='store_true', help='Enable debug output')
parser.add_argument('--no-progress', action='store_true', help='Disable progress bar during upload')
parser.add_argument('--no-reset', action='store_true', help='Skip device reset before uploading protocol (default: reset is enabled)')
parser.add_argument('--dry-run', action='store_true', help='Prepare protocol and show settings, but do not upload to device')
parser.add_argument('--no-calibration', action='store_true', help='Disable calibration matrix correction (calibration is enabled by default)')
parser.add_argument('--disable-current-scaling', action='store_true', help='Disable automatic current scaling - uses fixed 17mA current (default: current scaling is enabled)')
parser.add_argument('--verbose', action='store_true', help='Enable verbose protocol preparation output')
args = parser.parse_args()

device_id = args.device_id
selected_group = args.led_group

# Start with DEFAULT_OPTICAL_POWER as the base (4, 8, 12)
# This contains all 4 channels
power_matrix = DEFAULT_OPTICAL_POWER.copy()

# Handle legacy single-color mode
if args.power_matrix_csv:
	power_matrix *= 0
	# Legacy mode: set only one colour
	
	color = args.color
	channel = color if color < 4 else color - 4
	group = 1 if color < 4 else 2

	# Override just this channel with CSV data
	power_matrix[channel] = load_power_matrix_from_csv(args.power_matrix_csv)

	# Make sure this matches the selected group
	if group != selected_group:
		print(f"Warning: Color {color} uses LEDG{group}, but LEDG{selected_group} is selected. Switching to LEDG{group}.")
		selected_group = group

	print(f"Legacy mode: Loaded CSV for Color {color} → Channel {channel}, LEDG{group}")

print(f"Selected LED group: LEDG{selected_group}")

# Load calibration matrix
if args.no_calibration:
	print("Calibration disabled by --no-calibration flag")
	calibration_matrix, calibration_applied = None, False
else:
	calibration_matrix, calibration_applied = load_device_calibration(device_id, verbose=True, debug=args.debug)

# Identify which channels have non-zero power
active_channels = [ch for ch in range(4) if np.any(power_matrix[ch] > 0)]
if not active_channels:
	print("Error: No channels have non-zero power. Nothing to upload.")
	sys.exit(1)

print(f"Active channels: {active_channels}")

# PREPARE PROTOCOL
print("Preparing protocol...")

try:
	protocol = prepare_protocol(
		matrix_data=power_matrix,  # 3D array (4, 8, 12)
		on_time=args.on_time,
		period=args.period,
		blue_power_nominal=get_reference_powers(selected_group),
		start_time_matrix=np.zeros((NROWS, NCOLS)),
		on_time_matrix=np.full((NROWS, NCOLS), args.on_time),
		calibration_matrix=calibration_matrix,
		enable_current_scaling=not args.disable_current_scaling,
		lpwm_value=args.lpwm_value,
		verbose=args.verbose,
		led_group=selected_group,
		print_summary=True,
		device_id=device_id,
		temperature=args.temperature,
		verify=not args.no_verify,
		show_progress=not args.no_progress,
		reset_device=not args.no_reset
	)
except ValueError as e:
	print(f"Error preparing protocol:\n  {e}")
	sys.exit(1)

# Exit if dry-run mode
if args.dry_run:
	print("\nDRY RUN - Protocol not uploaded to device")
	sys.exit(0)

# UPLOAD PROTOCOL
# Create GUI and DIYA device
diya_device = diya(gui=None, device_id=device_id, debug=args.debug)
diya_device.start()

# Upload protocol
try:
	upload_protocol(diya_device, protocol)
except (ValueError, Exception) as e:
	if args.debug:
		import traceback
		traceback.print_exc()
	sys.exit(1)


