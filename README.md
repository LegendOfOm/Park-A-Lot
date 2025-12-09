# Park-A-Lot
SFML Environment + Algorithm
minGW must be in 64bit.
Install 64bit minGW using MSYS2 and then install sfml using MSYS2
Install SFML using MSYS2 (I insist you use MSYS instead of downloading it standalone since htis is much more convenient)
Exact Instructions:
1. Download MSYS2-x86_64: https://www.msys2.org/ 
2. Open MSYS2 and run: pacman -Syu
3. In MSYS2 (restart if necessary), Run: pacman -S mingw-w64-x86_64-gcc        pacman -S mingw-w64-x86_64-gdb
4. Now you should have C:\msys64\mingw64\bin. This needs to be added to path in environment variables.
5. You should now also have MSYS2 MINGW64 terminal as an accessible app. Open this terminal
6. Just like before, run pacman -Syu (restart after if necessary) and then pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-gdb
7. Run (In the same MSYS2 MINGW64 terminal) pacman -S mingw-w64-x86_64-sfml
8. Go through the .vscode Windows folder and rename it .vscode (make sure it is under SFML). Adjust all of the file paths IF NECESSARY (shouldn't be) in the .json to match your workspace.
9. That should be it. I already adjusted all of the vscode json files to work so you shouldn't have to go through the same suffering I did. Please open the sfml folder as the root path in the editor so that you don't have other .jsons messing with the setup.

MAC Instructions
1. Download Homebrew by opening terminal and run: /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
2. Add Homebrew to your shell environement. Normally this means these commands: echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile   eval "$(/opt/homebrew/bin/brew shellenv)"
3. Run in terminal: brew update    brew upgrade
4. Install the compiler (terminal. You probably already have this though): xcode-select --install
5. Instal sfml (terminal): brew install sfml
6. Adjust your IF NECESSARY (shouldn't be) in vsCode to match the Mac ones specifically.
7. For your folder layout, make sure you utilize the macOS main file nad macOS vscode. Rename the macOS main to main.cpp
