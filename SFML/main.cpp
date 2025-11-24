// main.cpp  --- SFML 3 particle playground
#include <SFML/Graphics.hpp>
#include <optional>
#include <random>
#include <vector>
#include <sstream>
#include <iomanip>
#include <iostream>

struct Particle
{
    sf::Vector2f position;
    sf::Vector2f velocity;
    float        radius;
    sf::Color    color;
    float        lifetime;   // seconds left
};

int main()
{
    // Window ------------------------------------------------------------------
    sf::RenderWindow window(sf::VideoMode({1280u, 720u}), "SFML 3 Playground");
    window.setVerticalSyncEnabled(true);

    // Font + FPS text ---------------------------------------------------------
    sf::Font font;
    if (!font.openFromFile("Ldfcomicsans-jj7l.ttf")) 
    {
        std::cerr << "Failed to load font\n";
        return 1;
    }

    sf::Text fpsText(font);
    fpsText.setCharacterSize(18);
    fpsText.setFillColor(sf::Color::White);
    fpsText.setPosition({10.f, 5.f});
    fpsText.setString("FPS: 0");

    // Random engine -----------------------------------------------------------
    std::mt19937 rng(static_cast<unsigned>(std::time(nullptr)));
    std::uniform_real_distribution<float> xDist(0.f, 1280.f);
    std::uniform_real_distribution<float> yDist(0.f, 720.f);
    std::uniform_real_distribution<float> velDist(-150.f, 150.f);
    std::uniform_real_distribution<float> radiusDist(4.f, 12.f);
    std::uniform_real_distribution<float> lifeDist(2.f, 6.f);

    auto randomColor = [&]() -> sf::Color {
        std::uniform_int_distribution<int> c(50, 255);
        return sf::Color(c(rng), c(rng), c(rng));
    };

    // Particle storage --------------------------------------------------------
    std::vector<Particle> particles;

    auto spawnParticleAt = [&](sf::Vector2f pos) {
        Particle p;
        p.position = pos;
        p.velocity = { velDist(rng), velDist(rng) - 50.f };
        p.radius   = radiusDist(rng);
        p.color    = randomColor();
        p.lifetime = lifeDist(rng);
        particles.push_back(p);
    };

    auto spawnBurst = [&](sf::Vector2f center, int count) {
        for (int i = 0; i < count; ++i)
            spawnParticleAt(center);
    };

    // Pre-spawn some particles
    for (int i = 0; i < 200; ++i)
        spawnParticleAt({xDist(rng), yDist(rng)});

    // Timing ------------------------------------------------------------------
    sf::Clock clock;
    float fpsSmoothing = 0.f;

    const float gravity = 400.f;

    // Main loop ---------------------------------------------------------------
    while (window.isOpen())
    {
        // --- Events ----------------------------------------------------------
        while (const std::optional event = window.pollEvent())
        {
            if (event->is<sf::Event::Closed>())
                window.close();

            if (const auto *key = event->getIf<sf::Event::KeyPressed>())
            {
                // Escape closes window
                if (key->scancode == sf::Keyboard::Scancode::Escape)
                    window.close();
                // Space spawns a burst in the center
                if (key->scancode == sf::Keyboard::Scancode::Space)
                    spawnBurst(sf::Vector2f{640.f, 360.f}, 80);
            }
        }

        // --- Delta time & FPS -----------------------------------------------
        float dt = clock.restart().asSeconds();
        // simple exponential smoothing for FPS readout
        float fps = dt > 0.f ? 1.f / dt : 0.f;
        fpsSmoothing = 0.9f * fpsSmoothing + 0.1f * fps;

        std::ostringstream ss;
        ss << "Particles: " << particles.size()
           << "   FPS: " << std::fixed << std::setprecision(1) << fpsSmoothing;
        fpsText.setString(ss.str());

        // --- Spawn with mouse (hold left button) ----------------------------
        if (sf::Mouse::isButtonPressed(sf::Mouse::Button::Left))
        {
            sf::Vector2i pixelPos = sf::Mouse::getPosition(window);
            sf::Vector2f worldPos{static_cast<float>(pixelPos.x),
                                  static_cast<float>(pixelPos.y)};
            spawnBurst(worldPos, 6);
        }

        // --- Update particles -----------------------------------------------
        const sf::Vector2f windowSize{
            static_cast<float>(window.getSize().x),
            static_cast<float>(window.getSize().y)
        };

        for (auto it = particles.begin(); it != particles.end(); )
        {
            Particle &p = *it;

            // physics
            p.velocity.y += gravity * dt;
            p.position   += p.velocity * dt;
            p.lifetime   -= dt;

            // bounce off walls
            if (p.position.x - p.radius < 0.f) {
                p.position.x = p.radius;
                p.velocity.x = std::abs(p.velocity.x);
            }
            if (p.position.x + p.radius > windowSize.x) {
                p.position.x = windowSize.x - p.radius;
                p.velocity.x = -std::abs(p.velocity.x);
            }
            if (p.position.y - p.radius < 0.f) {
                p.position.y = p.radius;
                p.velocity.y = std::abs(p.velocity.y);
            }
            if (p.position.y + p.radius > windowSize.y) {
                p.position.y = windowSize.y - p.radius;
                p.velocity.y = -std::abs(p.velocity.y) * 0.6f; // lose energy
            }

            // fade by lifetime
            float t = std::max(p.lifetime, 0.f);
            float alpha = std::min(t / 2.f, 1.f);  // quick fade near the end
            sf::Color c = p.color;
            c.a = static_cast<std::uint8_t>(alpha * 255.f);
            p.color = c;

            // remove dead
            if (p.lifetime <= 0.f || p.color.a == 0)
                it = particles.erase(it);
            else
                ++it;
        }

        // --- Render ---------------------------------------------------------
        window.clear(sf::Color(15, 15, 30));

        sf::CircleShape circle;
        for (const Particle &p : particles)
        {
            circle.setRadius(p.radius);
            circle.setOrigin({p.radius, p.radius});
            circle.setPosition(p.position);
            circle.setFillColor(p.color);
            window.draw(circle);
        }

        window.draw(fpsText);

        window.display();
    }

    return 0;
}
