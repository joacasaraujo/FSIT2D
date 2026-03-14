#pragma once
#include <iostream>
#include <cstdlib>

// Assert sempre ativo, inclusive em modo Release
#define ASSERT(expr) \
    do { \
        if (!(expr)) { \
            std::cerr << "❌ ASSERT FAILED: " << #expr << "\n" \
                      << "  → File: " << __FILE__ << "\n" \
                      << "  → Line: " << __LINE__ << "\n"; \
            std::abort(); \
        } \
    } while (0)

// Assert com mensagem explicativa
#define ASSERT_MSG(expr, msg) \
    do { \
        if (!(expr)) { \
            std::cerr << "❌ ASSERT FAILED: " << #expr << "\n" \
                      << "  → " << msg << "\n" \
                      << "  → File: " << __FILE__ << "\n" \
                      << "  → Line: " << __LINE__ << "\n"; \
            std::abort(); \
        } \
    } while (0)
    