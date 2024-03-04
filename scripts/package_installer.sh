#!/usr/bin/env bash

## Created: March 1, 2024
## Updated: March 1, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Checks installation of packages

PKG=$1

echo "Installing newest version of $PKG."

if command -v apt &> /dev/null; then
	sudo apt install -y "$PKG"
elif command -v yum &> /dev/null; then
	sudo yum install -y "$PKG"
elif command -v dnf &> /dev/null; then
	sudo dnf install -y "$PKG"
elif command -v zypper &> /dev/null; then
	sudo zypper --non-interactive install "$PKG"
elif command -v pacman &> /dev/null; then
	sudo pacman install --noconfirm "$PKG"
elif command -v brew &> /dev/null; then
	brew install -y "$PKG"
else
	echo "Can't determine package manager. Install $PKG manually."
	exit 1
fi
