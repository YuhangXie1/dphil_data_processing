"""
PNG IMAGE SEQUENCE -> LOSSLESS VIDEO
WITH INTERACTIVE MOUSE ROI SELECTION

Expected filenames:
    LED450NM_FOV4_frame_0000.png
    LED450NM_FOV4_frame_0001.png
    ...
    LED450NM_FOV4_frame_0125.png

FEATURES
--------
1. Select the start and end frame in the USER SETTINGS below.
2. Exclude individual frames.
3. Exclude ranges of frames.
4. Interactively select a rectangular ROI using the mouse.
5. Preview the selected ROI before making the video.
6. Keep the original pixel resolution of the selected ROI.
7. Create a lossless FFV1 video.
8. No resizing or lossy video compression is performed.

REQUIREMENTS
------------
Python packages:

    pip install opencv-python

FFmpeg must also be installed and available from the command line.

Windows:
    https://ffmpeg.org/download.html

macOS:
    brew install ffmpeg

Ubuntu/Debian:
    sudo apt install ffmpeg
"""

from pathlib import Path
import re
import subprocess
import tempfile

import cv2


# ============================================================
# USER SETTINGS
# ============================================================

# Folder containing the PNG images.
IMAGE_FOLDER = Path(
    r"FOV_no4_comparable_bright"
)

# Output video.
#
# MKV is recommended for FFV1 lossless video.
OUTPUT_VIDEO = Path(
    r"FOV_no4/LED450NM_FOV4_video_comparable_bright.mkv"
)


# ------------------------------------------------------------
# FILE NAME
# ------------------------------------------------------------

# Expected format:
#
# LED450NM_FOV4_frame_0000.png
# LED450NM_FOV4_frame_0001.png
# ...
#
FILE_PREFIX = "LED450NM_FOV4_frame_"


# ------------------------------------------------------------
# FRAME RANGE
# ------------------------------------------------------------

# Both values are INCLUSIVE.
#
# Example:
#
# START_FRAME = 20
# END_FRAME = 100
#
# means frames 0020 through 0100.

START_FRAME = 0
END_FRAME = 135


# ------------------------------------------------------------
# VIDEO FRAME RATE
# ------------------------------------------------------------

FPS = 5


# ============================================================
# EXCLUDE FRAMES
# ============================================================

# Individual frame numbers that should NOT appear in the video.
#
# Example:
#
# EXCLUDE_FRAMES = {
#     3,
#     7,
#     42,
# }

EXCLUDE_FRAMES = {
    18,
    19,
    24,
    30,
    31,
    33,
    34,
    42,
    50,
    53,
    54,
    58,
    60,
    79,
    80,
    81,
    82,
    83,
    84,
    85,
    86,
    87,
    92,
    93,
    108,
    113,
    114,
    118,
    119,
    120,
    126,
}


# ------------------------------------------------------------
# EXCLUDE RANGES
# ------------------------------------------------------------

# You can also exclude entire ranges.
#
# Each tuple is:
#
#     (first_frame, last_frame)
#
# Both numbers are inclusive.
#
# Example:
#
# EXCLUDE_RANGES = [
#     (20, 25),
#     (70, 75),
# ]
#
# removes:
#
# 20, 21, 22, 23, 24, 25
#
# and:
#
# 70, 71, 72, 73, 74, 75

EXCLUDE_RANGES = [
    # (20, 25),
    # (70, 75),
]


# ============================================================
# ROI SETTINGS
# ============================================================

# True:
#     A window opens and lets you draw the ROI with the mouse.
#
# False:
#     The complete image is used.
#
USE_INTERACTIVE_ROI = True


# Which frame should be displayed when selecting the ROI?
#
# None:
#     Use the first frame that will appear in the video.
#
# Or specify a frame number:
#
# ROI_PREVIEW_FRAME = 50
#
ROI_PREVIEW_FRAME = None


# ============================================================
# FUNCTIONS
# ============================================================


def get_frame_number(path):
    """
    Extract the numerical frame number from a PNG filename.

    Example:

        LED450NM_FOV4_frame_0042.png

    returns:

        42

    Files that do not match FILE_PREFIX are ignored.
    """

    pattern = rf"^{re.escape(FILE_PREFIX)}(\d+)\.png$"

    match = re.match(
        pattern,
        path.name,
        flags=re.IGNORECASE
    )

    if match is None:
        return None

    return int(match.group(1))


def build_exclusion_set():
    """
    Combine individual excluded frames and excluded ranges
    into one Python set.
    """

    excluded = set(EXCLUDE_FRAMES)

    for start, end in EXCLUDE_RANGES:

        if end < start:
            raise ValueError(
                f"Invalid exclusion range: ({start}, {end})\n"
                "The end must be >= the start."
            )

        excluded.update(
            range(start, end + 1)
        )

    return excluded


def find_selected_images():
    """
    Search IMAGE_FOLDER and determine exactly which images
    should appear in the video.

    Returns:

        [
            (frame_number, Path),
            (frame_number, Path),
            ...
        ]

    sorted numerically by frame number.
    """

    if not IMAGE_FOLDER.exists():

        raise FileNotFoundError(
            f"Image folder does not exist:\n"
            f"{IMAGE_FOLDER}"
        )

    if END_FRAME < START_FRAME:

        raise ValueError(
            "END_FRAME must be >= START_FRAME."
        )

    excluded = build_exclusion_set()

    images = []

    for path in IMAGE_FOLDER.iterdir():

        # Ignore folders.
        if not path.is_file():
            continue

        frame_number = get_frame_number(path)

        # Ignore unrelated files.
        if frame_number is None:
            continue

        # Ignore images outside the requested range.
        if not (
            START_FRAME
            <= frame_number
            <= END_FRAME
        ):
            continue

        # Ignore explicitly excluded frames.
        if frame_number in excluded:
            continue

        images.append(
            (frame_number, path)
        )

    # Sort by actual frame number.
    images.sort(
        key=lambda item: item[0]
    )

    return images


def read_image(path):
    """
    Read a PNG while preserving its original image format
    as much as OpenCV allows.

    IMREAD_UNCHANGED prevents OpenCV from automatically
    forcing the image to standard 8-bit RGB.
    """

    image = cv2.imread(
        str(path),
        cv2.IMREAD_UNCHANGED
    )

    if image is None:

        raise RuntimeError(
            f"Could not read image:\n{path}"
        )

    return image


def make_roi_display_image(image):
    """
    Create an image suitable for OpenCV's ROI selector.

    IMPORTANT:

    This image is ONLY used for the interactive preview.

    The actual video is generated directly from the original
    PNG files using FFmpeg.

    Therefore, changing the preview representation here does
    NOT reduce the quality of the final video.
    """

    # --------------------------------------------------------
    # 16-bit images
    # --------------------------------------------------------
    #
    # OpenCV's display window generally needs an 8-bit image
    # for convenient visualization.
    #
    # If the original PNG is 16-bit, scale it to 8-bit ONLY
    # for display.
    #
    # The original 16-bit file remains untouched.
    # --------------------------------------------------------

    if image.dtype != "uint8":

        minimum = image.min()
        maximum = image.max()

        if maximum > minimum:

            display = (
                (image.astype("float32") - minimum)
                / (maximum - minimum)
                * 255
            ).astype("uint8")

        else:

            display = image.astype("uint8")

    else:

        display = image.copy()

    return display


def select_roi_interactively(image):
    """
    Open an interactive ROI-selection window.

    If the source image is larger than the desired preview size,
    the DISPLAY COPY is scaled down so the entire image fits on
    screen.

    IMPORTANT:
    The original image is NOT resized.

    After the user selects an ROI on the scaled preview, the
    coordinates are converted back to the ORIGINAL full-resolution
    image coordinates.

    Therefore the final video still uses the original image pixels.

    Controls
    --------
    LEFT MOUSE:
        Click and drag to select the ROI.

    ENTER or SPACE:
        Accept the selected ROI.

    C:
        Cancel.
    """

    # ========================================================
    # PREVIEW WINDOW SIZE
    # ========================================================
    #
    # Maximum size of the ROI selection image on your screen.
    #
    # Change these if desired.
    #
    # For example, on a 1920 x 1080 monitor:
    #
    #     1400 x 850
    #
    # leaves room for the window title bar, taskbar, etc.
    #
    # This ONLY affects the preview.
    # It does NOT affect the final video resolution.

    MAX_PREVIEW_WIDTH = 1400
    MAX_PREVIEW_HEIGHT = 850

    # --------------------------------------------------------
    # Prepare image for display
    # --------------------------------------------------------

    display = make_roi_display_image(image)

    original_height, original_width = image.shape[:2]

    print("\n" + "=" * 60)
    print("INTERACTIVE ROI SELECTION")
    print("=" * 60)

    print(
        f"\nOriginal image resolution: "
        f"{original_width} x {original_height}"
    )

    # --------------------------------------------------------
    # Determine preview scaling
    # --------------------------------------------------------
    #
    # We calculate how much the image needs to be scaled to fit
    # inside MAX_PREVIEW_WIDTH x MAX_PREVIEW_HEIGHT.
    #
    # scale = 1.0 means no scaling.
    #
    # scale < 1.0 means the preview is reduced.
    #
    # We NEVER scale above 1.0 because there is no reason to
    # enlarge small images for ROI selection.
    # --------------------------------------------------------

    width_scale = (
        MAX_PREVIEW_WIDTH / original_width
    )

    height_scale = (
        MAX_PREVIEW_HEIGHT / original_height
    )

    scale = min(
        width_scale,
        height_scale,
        1.0
    )

    # --------------------------------------------------------
    # Resize DISPLAY COPY if necessary
    # --------------------------------------------------------

    if scale < 1.0:

        preview_width = round(
            original_width * scale
        )

        preview_height = round(
            original_height * scale
        )

        display_scaled = cv2.resize(
            display,
            (preview_width, preview_height),
            interpolation=cv2.INTER_AREA
        )

        print(
            f"Preview resolution: "
            f"{preview_width} x {preview_height}"
        )

        print(
            f"Preview scale: "
            f"{scale:.4f}"
        )

    else:

        display_scaled = display

        preview_height, preview_width = (
            display_scaled.shape[:2]
        )

        print(
            "\nImage already fits inside the preview window."
        )

    # --------------------------------------------------------
    # Instructions
    # --------------------------------------------------------

    print(
        "\nDraw a rectangle using the LEFT MOUSE BUTTON."
    )

    print(
        "Press ENTER or SPACE to accept."
    )

    print(
        "Press C to cancel."
    )

    print(
        "\nNOTE: The displayed image may be scaled down, "
        "but the final ROI is calculated using the ORIGINAL "
        "full-resolution image."
    )

    # --------------------------------------------------------
    # Interactive selection
    # --------------------------------------------------------

    roi_scaled = cv2.selectROI(
        "Select ROI - ENTER/SPACE = Accept, C = Cancel",
        display_scaled,
        showCrosshair=True,
        fromCenter=False
    )

    cv2.destroyAllWindows()

    # ROI returned by OpenCV refers to the SCALED preview.
    scaled_x, scaled_y, scaled_w, scaled_h = [
        int(value)
        for value in roi_scaled
    ]

    # --------------------------------------------------------
    # Detect cancellation
    # --------------------------------------------------------

    if scaled_w == 0 or scaled_h == 0:

        raise RuntimeError(
            "\nROI selection was cancelled or the selected "
            "area had zero width/height."
        )

    # ========================================================
    # CONVERT BACK TO ORIGINAL IMAGE COORDINATES
    # ========================================================
    #
    # Example:
    #
    # Original:
    #
    #     4000 x 3000
    #
    # Preview:
    #
    #     1000 x 750
    #
    # scale = 0.25
    #
    # If you select:
    #
    #     x = 100
    #
    # on the preview, the corresponding original coordinate is:
    #
    #     100 / 0.25 = 400
    #
    # Therefore we divide all coordinates by the scale.
    # --------------------------------------------------------

    x = round(
        scaled_x / scale
    )

    y = round(
        scaled_y / scale
    )

    w = round(
        scaled_w / scale
    )

    h = round(
        scaled_h / scale
    )

    # --------------------------------------------------------
    # Keep coordinates inside the original image.
    # --------------------------------------------------------

    x = max(
        0,
        min(x, original_width - 1)
    )

    y = max(
        0,
        min(y, original_height - 1)
    )

    # Make sure width/height don't extend outside the image.

    w = min(
        w,
        original_width - x
    )

    h = min(
        h,
        original_height - y
    )

    # --------------------------------------------------------
    # Final validation
    # --------------------------------------------------------

    if w <= 0 or h <= 0:

        raise RuntimeError(
            "The calculated ROI has zero width or height."
        )

    # --------------------------------------------------------
    # Print selection information
    # --------------------------------------------------------

    print("\n" + "=" * 60)
    print("ROI CONVERTED TO ORIGINAL IMAGE COORDINATES")
    print("=" * 60)

    print(
        f"\nPreview selection:"
        f"\n  X      = {scaled_x}"
        f"\n  Y      = {scaled_y}"
        f"\n  Width  = {scaled_w}"
        f"\n  Height = {scaled_h}"
    )

    print(
        f"\nOriginal-resolution selection:"
        f"\n  X      = {x}"
        f"\n  Y      = {y}"
        f"\n  Width  = {w}"
        f"\n  Height = {h}"
    )

    print(
        f"\nFinal video resolution:"
        f"\n  {w} x {h} pixels"
    )

    return x, y, w, h


def show_selected_roi_preview(image, roi):
    """
    Show ONLY the selected ROI at its actual pixel resolution.

    This gives you a chance to visually inspect the area
    before the video is generated.

    Press any key to continue.
    """

    x, y, w, h = roi

    cropped = image[
        y:y + h,
        x:x + w
    ]

    display = make_roi_display_image(
        cropped
    )

    print("\nROI PREVIEW")
    print("-" * 60)

    print(f"X      : {x}")
    print(f"Y      : {y}")
    print(f"Width  : {w}")
    print(f"Height : {h}")

    print(
        f"\nFinal video resolution will be:"
        f"\n{w} x {h} pixels"
    )

    print(
        "\nA preview of the selected area will now open."
    )

    print(
        "Press any key while the preview window is active "
        "to continue."
    )

    cv2.imshow(
        "Selected ROI Preview - Press Any Key To Continue",
        display
    )

    cv2.waitKey(0)

    cv2.destroyAllWindows()


def find_roi_preview_image(images):
    """
    Determine which image should be displayed for selecting
    the ROI.

    If ROI_PREVIEW_FRAME is None:
        use the first selected video frame.

    Otherwise:
        attempt to find the specified frame.
    """

    if ROI_PREVIEW_FRAME is None:

        return images[0]

    # Search all PNG files because the requested preview frame
    # does not necessarily have to be part of the final video.

    for path in IMAGE_FOLDER.iterdir():

        if not path.is_file():
            continue

        frame_number = get_frame_number(path)

        if frame_number == ROI_PREVIEW_FRAME:

            return frame_number, path

    raise RuntimeError(
        f"\nROI preview frame "
        f"{ROI_PREVIEW_FRAME:04d} "
        f"could not be found."
    )


def validate_images(images):
    """
    Verify that all selected images have the same:

        - width
        - height
        - number of channels
        - data type / bit depth

    This prevents unexpected changes during video generation.
    """

    first_number, first_path = images[0]

    first_image = read_image(
        first_path
    )

    reference_shape = first_image.shape
    reference_dtype = first_image.dtype

    print("\nChecking all selected images...")

    for index, (
        frame_number,
        path
    ) in enumerate(
        images,
        start=1
    ):

        image = read_image(path)

        # Check dimensions/channels.

        if image.shape != reference_shape:

            raise ValueError(
                f"\n\nImage shape mismatch at "
                f"frame {frame_number:04d}:\n\n"
                f"{path}\n\n"
                f"Expected: {reference_shape}\n"
                f"Found:    {image.shape}"
            )

        # Check image data type.

        if image.dtype != reference_dtype:

            raise ValueError(
                f"\n\nImage data type mismatch at "
                f"frame {frame_number:04d}:\n\n"
                f"{path}\n\n"
                f"Expected: {reference_dtype}\n"
                f"Found:    {image.dtype}"
            )

        print(
            f"\rChecking image "
            f"{index}/{len(images)}",
            end="",
            flush=True
        )

    print("\nImage check complete.")

    return first_image


def create_lossless_video(
    images,
    roi=None
):
    """
    Create the video using FFmpeg and the FFV1 lossless codec.

    Parameters
    ----------
    images:
        List of selected image paths.

    roi:
        None

        OR:

        (x, y, width, height)

    The PNG images are fed directly to FFmpeg.

    If an ROI is selected, FFmpeg performs the crop directly.

    There is:

        - NO resizing
        - NO JPEG conversion
        - NO H.264 compression
        - NO H.265 compression

    FFV1 is mathematically lossless.
    """

    with tempfile.TemporaryDirectory() as temp_dir:

        temp_dir = Path(temp_dir)

        concat_file = (
            temp_dir / "frames.txt"
        )

        # ----------------------------------------------------
        # Build FFmpeg's list of source PNG files.
        # ----------------------------------------------------

        with open(
            concat_file,
            "w",
            encoding="utf-8"
        ) as f:

            for frame_number, path in images:

                absolute_path = (
                    path.resolve()
                )

                # Escape single quotes for FFmpeg's concat
                # demuxer syntax.

                safe_path = str(
                    absolute_path
                ).replace(
                    "'",
                    "'\\''"
                )

                f.write(
                    f"file '{safe_path}'\n"
                )

                f.write(
                    f"duration "
                    f"{1.0 / FPS:.12f}\n"
                )

            # Repeat final frame so FFmpeg honors the final
            # frame duration.

            last_path = (
                images[-1][1].resolve()
            )

            safe_last_path = str(
                last_path
            ).replace(
                "'",
                "'\\''"
            )

            f.write(
                f"file '{safe_last_path}'\n"
            )

        # ----------------------------------------------------
        # Base FFmpeg command
        # ----------------------------------------------------

        command = [
            "ffmpeg",

            # Automatically overwrite an existing output.
            "-y",

            # Input is an FFmpeg concat list.
            "-f",
            "concat",

            # Allow absolute filenames.
            "-safe",
            "0",

            # Input list.
            "-i",
            str(concat_file),
        ]

        # ----------------------------------------------------
        # ROI CROP
        # ----------------------------------------------------

        if roi is not None:

            x, y, w, h = roi

            crop_filter = (
                f"crop="
                f"{w}:"
                f"{h}:"
                f"{x}:"
                f"{y}"
            )

            command.extend(
                [
                    "-vf",
                    crop_filter
                ]
            )

        # ----------------------------------------------------
        # LOSSLESS VIDEO SETTINGS
        # ----------------------------------------------------

        command.extend(
            [
                # Desired playback rate.
                "-r",
                str(FPS),

                # FFV1 lossless video codec.
                "-c:v",
                "ffv1",

                # FFV1 version/level 3.
                "-level",
                "3",

                # Add CRC information for slices.
                "-slicecrc",
                "1",

                # Matroska container.
                "-f",
                "matroska",

                str(OUTPUT_VIDEO)
            ]
        )

        print("\n" + "=" * 60)
        print("CREATING VIDEO")
        print("=" * 60)

        print("\nFFmpeg command:\n")

        print(
            " ".join(
                f'"{part}"'
                for part in command
            )
        )

        print()

        try:

            subprocess.run(
                command,
                check=True
            )

        except FileNotFoundError:

            raise RuntimeError(
                "\nFFmpeg could not be found.\n\n"
                "Make sure FFmpeg is installed and that "
                "the 'ffmpeg' command is available from "
                "your command line."
            )

        except subprocess.CalledProcessError as error:

            raise RuntimeError(
                "\nFFmpeg failed while creating the video."
            ) from error


# ============================================================
# MAIN PROGRAM
# ============================================================


def main():

    print("=" * 60)
    print("PNG -> LOSSLESS VIDEO")
    print("=" * 60)

    # --------------------------------------------------------
    # Find requested images.
    # --------------------------------------------------------

    images = find_selected_images()

    if not images:

        raise RuntimeError(
            "\nNo images were found matching the current "
            "settings.\n\n"
            "Check:\n"
            "  IMAGE_FOLDER\n"
            "  FILE_PREFIX\n"
            "  START_FRAME\n"
            "  END_FRAME\n"
            "  EXCLUDE_FRAMES\n"
            "  EXCLUDE_RANGES"
        )

    print("\nImage folder:")
    print(IMAGE_FOLDER)

    print("\nRequested frame range:")

    print(
        f"{START_FRAME:04d} -> "
        f"{END_FRAME:04d}"
    )

    print(
        f"\nNumber of frames going into video: "
        f"{len(images)}"
    )

    print("\nFirst selected frame:")

    print(
        f"{images[0][0]:04d} : "
        f"{images[0][1].name}"
    )

    print("\nLast selected frame:")

    print(
        f"{images[-1][0]:04d} : "
        f"{images[-1][1].name}"
    )

    # --------------------------------------------------------
    # Report excluded frames.
    # --------------------------------------------------------

    excluded = build_exclusion_set()

    excluded_in_range = sorted(
        frame
        for frame in excluded
        if (
            START_FRAME
            <= frame
            <= END_FRAME
        )
    )

    if excluded_in_range:

        print(
            "\nExplicitly excluded frames:"
        )

        print(excluded_in_range)

    else:

        print(
            "\nNo frames explicitly excluded."
        )

    # --------------------------------------------------------
    # Validate image sequence.
    # --------------------------------------------------------

    first_image = validate_images(
        images
    )

    height, width = (
        first_image.shape[:2]
    )

    print("\nSource image properties:")

    print(
        f"Resolution : "
        f"{width} x {height}"
    )

    print(
        f"Data type  : "
        f"{first_image.dtype}"
    )

    if first_image.ndim == 2:

        print(
            "Channels   : 1"
        )

    else:

        print(
            f"Channels   : "
            f"{first_image.shape[2]}"
        )

    # --------------------------------------------------------
    # INTERACTIVE ROI SELECTION
    # --------------------------------------------------------

    roi = None

    if USE_INTERACTIVE_ROI:

        preview_frame_number, preview_path = (
            find_roi_preview_image(
                images
            )
        )

        print(
            f"\nUsing frame "
            f"{preview_frame_number:04d} "
            f"for ROI selection."
        )

        print(
            preview_path
        )

        preview_image = read_image(
            preview_path
        )

        roi = select_roi_interactively(
            preview_image
        )

        # Print selected coordinates.

        x, y, w, h = roi

        print("\n" + "=" * 60)
        print("SELECTED ROI")
        print("=" * 60)

        print(
            f"\nTop-left X : {x}"
        )

        print(
            f"Top-left Y : {y}"
        )

        print(
            f"Width      : {w}"
        )

        print(
            f"Height     : {h}"
        )

        print(
            f"\nBottom-right coordinate:"
        )

        print(
            f"X = {x + w}"
        )

        print(
            f"Y = {y + h}"
        )

        print(
            f"\nVideo resolution:"
        )

        print(
            f"{w} x {h} pixels"
        )

        # ----------------------------------------------------
        # Show cropped result before making video.
        # ----------------------------------------------------

        show_selected_roi_preview(
            preview_image,
            roi
        )

    else:

        print(
            "\nInteractive ROI disabled."
        )

        print(
            "The complete image will be used."
        )

        print(
            f"\nVideo resolution:"
            f"\n{width} x {height}"
        )

    # --------------------------------------------------------
    # Make output directory if necessary.
    # --------------------------------------------------------

    OUTPUT_VIDEO.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    # --------------------------------------------------------
    # CREATE VIDEO
    # --------------------------------------------------------

    create_lossless_video(
        images,
        roi
    )

    # --------------------------------------------------------
    # DONE
    # --------------------------------------------------------

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)

    print(
        "\nVideo written to:"
    )

    print(
        OUTPUT_VIDEO
    )

    print(
        f"\nNumber of video frames: "
        f"{len(images)}"
    )

    print(
        f"Frame rate: "
        f"{FPS} FPS"
    )

    duration = (
        len(images) / FPS
    )

    print(
        f"Approximate duration: "
        f"{duration:.3f} seconds"
    )

    if roi is not None:

        x, y, w, h = roi

        print(
            f"\nROI:"
            f"\n  X      = {x}"
            f"\n  Y      = {y}"
            f"\n  Width  = {w}"
            f"\n  Height = {h}"
        )

        print(
            f"\nFinal resolution: "
            f"{w} x {h}"
        )

    else:

        print(
            f"\nFinal resolution: "
            f"{width} x {height}"
        )


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":
    main()