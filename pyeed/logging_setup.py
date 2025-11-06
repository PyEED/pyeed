import sys

from loguru import logger


def setup_logging() -> None:
    logger.remove()

    logger.add(
        sink=sys.stderr,
        level="WARNING",
        backtrace=True,
        diagnose=False,
        format="<green>{time:HH:mm:ss.SSS}</green> | "
        "<level>{level: <8}</level> | "
        "{message} | {extra}",
    )

    logger.add(
        "pyeed.log",
        rotation="10 MB",
        retention=5,
        compression="gz",
        level="DEBUG",
        format=(
            "{time:YYYY-MM-DD HH:mm:ss.SSS} | {level: <8} | {name}:{function}:{line} - {message}"
        ),
    )
