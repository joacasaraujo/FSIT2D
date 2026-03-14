#pragma once
#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>

class Timer {
public:
    // rótulo opcional (aparece na impressão)
    explicit Timer(std::string label = "Simulation Time")
        : _label(std::move(label)),
          _start(std::chrono::steady_clock::now()) {}

    // garante que não haja dupla impressão
    ~Timer() { if (!_stopped) Stop(); }

    // reinicia contagem sem imprimir
    void reset(std::string new_label = "") {
        if (!new_label.empty()) _label = std::move(new_label);
        _start = std::chrono::steady_clock::now();
        _stopped = false;
    }

    // imprime no formato hh:mm:ss.ms
    void Stop() {
        if (_stopped) return;
        _stopped = true;

        const auto end = std::chrono::steady_clock::now();
        const auto ms  = std::chrono::duration_cast<std::chrono::milliseconds>(end - _start).count();

        const auto h   = ms / 3600000;
        const auto m   = (ms % 3600000) / 60000;
        const auto s   = (ms % 60000) / 1000;
        const auto rem = ms % 1000;

        std::cout << "------------------------------------------\n";
        std::cout << _label << ": "
                  << std::setfill('0') << std::setw(2) << h << ":"
                  << std::setw(2) << m << ":"
                  << std::setw(2) << s << "."
                  << std::setw(3) << rem << " (hh:mm:ss.ms)\n";
    }

    // tempo decorrido (ms) sem encerrar
    long long elapsed_ms() const {
        const auto now = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(now - _start).count();
    }

private:
    std::string _label;
    std::chrono::time_point<std::chrono::steady_clock> _start;
    bool _stopped = false;
};
