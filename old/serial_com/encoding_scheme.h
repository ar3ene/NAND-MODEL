/**************************************************************************
 * @copyright 2020 EIGENBIT. All Rights Reserverd.
 *
 * @file host_transmitter.h
 *
 * @brief defines the structure for the host transmitter.
 *
 * @author(s) Wang Yingke,
 **************************************************************************/

#pragma once

#ifndef _ENCODING_SCHEME_H_
#define _ENCODING_SCHEME_H_

#include <windows.h>
#include <vector>

void encode_data(std::vector<BYTE>& data_in, std::vector<BYTE>& encoded_data,
                 int original_code_size, int target_code_size);

void decode_data(std::vector<BYTE>& encoded_data, std::vector<BYTE>& decoded_data,
                 int original_code_size, int target_code_size);

#endif // !_ENCODING_SCHEME_H_
