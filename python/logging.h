#ifndef PYMAXFLOW_PYTHON_LOGGING_H__
#define PYMAXFLOW_PYTHON_LOGGING_H__

#include <util/Logger.h>

namespace pymaxflow {

extern logger::LogChannel pylog;

/**
 * Get the log level of the python wrappers.
 */
logger::LogLevel getLogLevel();

/**
 * Set the log level of the python wrappers.
 */
void setLogLevel(logger::LogLevel logLevel);

} // namespace pycmc

#endif // PYMAXFLOW_PYTHON_LOGGING_H__

