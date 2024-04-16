# Fluidity Calculator for ImageJ/Fiji

This is an ImageJ/Fiji macro tool to calculate membrane fluidity values (aka GP) from spectral imaging or two-channel microscopy data.
Developed by [Pablo Carravilla](https://orcid.org/0000-0001-6592-7630) while working at the [CSI:Nano lab](https://www.csi-nano.org) at SciLifeLab and the Karolinska Institute.

The tool is a user-friendly graphical interface to calculate GP values directly in ImageJ/Fiji.
A Fiji installation is required, Fiji can be downloaded [here](https://fiji.sc).

The tool is installed by copying the files in the Fiji installation directory.
After installation, the tool is loaded from the 'More Tools' menu (double arrow icon, >>).
It can then be directly executed directly by clicking the icon on the toolbar.

## Installation Instructions
1. Find your Fiji installation directory and go to the _macros/toolsets_ folder.
	Mac OS:  
		Right click on Fiji in the Applications folder -> Show Package Contents
		or
		Go to this path: _/Applications/Fiji.app/macros/toolsets_
	Windows:
		It depends where you extracted the program when you downloaded it.

2. Delete any older version if you installed one in the past.

3. Copy the macro file 'FluidityCalculator_1.0.ijm' in the _toolsets_ folder.

4. Copy the 'icons' folder in the _toolsets_ folder.

5. (Optional) Copy the 'gp-viridis.lut' in the Fiji directory _luts_ folder.
	Mac OS: /Applications/Fiji.app/luts

6. Open Fiji.

7. Click in the More Tools button on the Toolbar, that is the right most icon of a double arrow (>>).

8. A new icon will appear in your toolbar, after opening an image click on it to launch the tool.
