#pragma once

//! Demange GCC-mangled typenames
//! \param name A mangled type-name
//! \return A demangled type-name
const char* demangle(const char* name);