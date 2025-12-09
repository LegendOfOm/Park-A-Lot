#pragma once
#include <SFML/Graphics.hpp>
#include <optional>
#include <string>

// Values the user chooses in the UI
struct UserConfig {
    double lotWidth;   // metres (x)
    double lotLength;  // metres (y)
    double stallWidth; // along lane
    double stallDepth; // perpendicular to lane
};

// Try a bunch of font locations so it works on Windows/macOS/Linux
bool loadSafeFont(sf::Font& font);

// Main UI entry point – opens a window, returns settings or std::nullopt if cancelled
std::optional<UserConfig> runUI();
