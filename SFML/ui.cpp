#include "ui.hpp"
#include <SFML/Window/Event.hpp>
#include <iostream>
#include <charconv>

// ---------------- Font loading (cross-platform) ----------------

bool loadSafeFont(sf::Font& font) {
    bool loaded = false;

#if defined(_WIN32)
    if (!loaded) loaded = font.openFromFile("C:\\Windows\\Fonts\\comic.ttf");
    if (!loaded) loaded = font.openFromFile("C:\\Windows\\Fonts\\comicbd.ttf");
#endif

    if (!loaded) loaded = font.openFromFile("assets/comic.ttf");

#if defined(__APPLE__)
    if (!loaded) loaded = font.openFromFile("/System/Library/Fonts/Supplemental/Arial.ttf");
    if (!loaded) loaded = font.openFromFile("/System/Library/Fonts/Helvetica.ttc");
#endif

#if defined(__linux__)
    if (!loaded) loaded = font.openFromFile("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf");
#endif

    if (!loaded) {
        std::cerr << "WARNING: could not load Comic Sans – trying assets/arial.ttf\n";
        font.openFromFile("assets/arial.ttf");
    }

    return true; 
}

// ---------------- Simple numeric input box ----------------

class InputBox {
public:
    sf::RectangleShape box;
    sf::Text text;
    std::string value;
    bool active = false;

    InputBox(float x, float y, const sf::Font& font, const std::string& placeholder = "")
        : box({200.f, 40.f})
        , text(font, placeholder, 20)
    {
        box.setPosition({x, y});
        box.setFillColor(sf::Color(230, 230, 230));
        box.setOutlineThickness(2.f);
        box.setOutlineColor(sf::Color(120, 120, 120));

        text.setFillColor(sf::Color::Black);
        text.setPosition({x + 5.f, y + 7.f});
    }

    void setActive(bool on) {
        active = on;
        box.setOutlineColor(on ? sf::Color::Blue : sf::Color(120, 120, 120));
    }

    void handleText(char32_t code) {
        if (!active) return;

        // Backspace
        if (code == U'\b' || code == 8) {
            if (!value.empty()) {
                value.pop_back();
                text.setString(value);
            }
            return;
        }

        // Digits or one '.' allowed
        if ((code >= U'0' && code <= U'9') || code == U'.') {
            if (code == U'.' && value.find('.') != std::string::npos) {
                return;
            }
            value += static_cast<char>(code);
            text.setString(value);
        }
    }

    bool clicked(float mx, float my) const {
        return box.getGlobalBounds().contains({mx, my});
    }

    std::optional<double> asDouble() const {
        if (value.empty()) return std::nullopt;

        double d;
        auto* b = value.data();
        auto* e = value.data() + value.size();
        if (std::from_chars(b, e, d).ec == std::errc()) return d;

        try { return std::stod(value); }
        catch (...) { return std::nullopt; }
    }

    void setValue(double v) {
        value = std::to_string(v);
        while (!value.empty() && value.back() == '0') value.pop_back();
        if (!value.empty() && value.back() == '.') value.pop_back();
        text.setString(value);
    }
};

// ---------------- Main UI window ----------------

std::optional<UserConfig> runUI() {
    sf::RenderWindow window(sf::VideoMode({700u, 430u}), "Parking Lot Setup");
    window.setFramerateLimit(60);

    sf::Font font;
    loadSafeFont(font);

    sf::Text title(font, "Parking Lot Simulator", 32);
    title.setFillColor(sf::Color::White);
    title.setPosition({190.f, 15.f});

    // Labels
    sf::Text lotLabel(font, "Lot dimensions (metres):", 20);
    lotLabel.setFillColor(sf::Color::White);
    lotLabel.setPosition({40.f, 80.f});

    sf::Text stallLabel(font, "Stall size (metres):", 20);
    stallLabel.setFillColor(sf::Color::White);
    stallLabel.setPosition({40.f, 180.f});

    sf::Text widthLabel(font, "Width (x)", 18);
    widthLabel.setFillColor(sf::Color::White);
    widthLabel.setPosition({40.f, 110.f});

    sf::Text lengthLabel(font, "Length (y)", 18);
    lengthLabel.setFillColor(sf::Color::White);
    lengthLabel.setPosition({40.f, 140.f});

    sf::Text stallWLabel(font, "Stall width (frontage)", 18);
    stallWLabel.setFillColor(sf::Color::White);
    stallWLabel.setPosition({40.f, 210.f});

    sf::Text stallDLabel(font, "Stall depth", 18);
    stallDLabel.setFillColor(sf::Color::White);
    stallDLabel.setPosition({40.f, 240.f});

    // Input boxes
    InputBox lotWidthBox (260.f, 105.f, font, "15");
    InputBox lotLengthBox(260.f, 135.f, font, "50");
    InputBox stallWBox   (260.f, 205.f, font, "2.5");
    InputBox stallDBox   (260.f, 235.f, font, "5.0");

    lotWidthBox.setValue(15.0);
    lotLengthBox.setValue(50.0);
    stallWBox.setValue(2.5);
    stallDBox.setValue(5.0);

    // Start button
    sf::RectangleShape startBtn({180.f, 50.f});
    startBtn.setPosition({40.f, 320.f});
    startBtn.setFillColor(sf::Color(120, 200, 120));

    sf::Text startText(font, "Generate Layout", 22);
    startText.setFillColor(sf::Color::Black);
    startText.setPosition({50.f, 330.f});

    // Initially focus lot width
    lotWidthBox.setActive(true);

    while (window.isOpen()) {

        while (const std::optional<sf::Event> ev = window.pollEvent()) {

            if (ev->is<sf::Event::Closed>()) {
                window.close();
                return std::nullopt;
            }

            if (const auto* te = ev->getIf<sf::Event::TextEntered>()) {
                if (lotWidthBox.active)  lotWidthBox.handleText(te->unicode);
                if (lotLengthBox.active) lotLengthBox.handleText(te->unicode);
                if (stallWBox.active)    stallWBox.handleText(te->unicode);
                if (stallDBox.active)    stallDBox.handleText(te->unicode);
            }

            if (const auto* mb = ev->getIf<sf::Event::MouseButtonPressed>()) {
                float mx = float(mb->position.x);
                float my = float(mb->position.y);

                // Clicking input boxes
                if (lotWidthBox.clicked(mx, my)) {
                    lotWidthBox.setActive(true);
                    lotLengthBox.setActive(false);
                    stallWBox.setActive(false);
                    stallDBox.setActive(false);
                } 
                else if (lotLengthBox.clicked(mx, my)) {
                    lotWidthBox.setActive(false);
                    lotLengthBox.setActive(true);
                    stallWBox.setActive(false);
                    stallDBox.setActive(false);
                } 
                else if (stallWBox.clicked(mx, my)) {
                    lotWidthBox.setActive(false);
                    lotLengthBox.setActive(false);
                    stallWBox.setActive(true);
                    stallDBox.setActive(false);
                } 
                else if (stallDBox.clicked(mx, my)) {
                    lotWidthBox.setActive(false);
                    lotLengthBox.setActive(false);
                    stallWBox.setActive(false);
                    stallDBox.setActive(true);
                } 
                else {
                    lotWidthBox.setActive(false);
                    lotLengthBox.setActive(false);
                    stallWBox.setActive(false);
                    stallDBox.setActive(false);
                }

                // Start button click
                if (startBtn.getGlobalBounds().contains({mx, my})) {
                    auto w  = lotWidthBox.asDouble();
                    auto l  = lotLengthBox.asDouble();
                    auto sw = stallWBox.asDouble();
                    auto sd = stallDBox.asDouble();

                    if (w && l && sw && sd &&
                        *w  > 0.0 && *l  > 0.0 &&
                        *sw > 0.0 && *sd > 0.0)
                    {
                        UserConfig cfg{*w, *l, *sw, *sd};
                        window.close();
                        return cfg;
                    } else {
                        std::cerr << "Invalid or empty values.\n";
                    }
                }
            }
        }

        // DRAW
        window.clear(sf::Color(40, 40, 40));

        window.draw(title);
        window.draw(lotLabel);
        window.draw(stallLabel);
        window.draw(widthLabel);
        window.draw(lengthLabel);
        window.draw(stallWLabel);
        window.draw(stallDLabel);

        window.draw(lotWidthBox.box);
        window.draw(lotWidthBox.text);
        window.draw(lotLengthBox.box);
        window.draw(lotLengthBox.text);
        window.draw(stallWBox.box);
        window.draw(stallWBox.text);
        window.draw(stallDBox.box);
        window.draw(stallDBox.text);

        window.draw(startBtn);
        window.draw(startText);

        window.display();
    }

    return std::nullopt;
}
