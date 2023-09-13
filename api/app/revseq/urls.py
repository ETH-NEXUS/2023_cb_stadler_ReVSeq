from django.contrib import admin
from django.urls import path, include, re_path
from users.views import CsrfCookieView, LoginView, LogoutView, UserViewSet
from core.views import (
    SampleCountViewSet,
    MetadataViewSet,
    PlateViewSet,
    SubstrainViewSet,
    SampleViewSet,
    download_file,
    FileViewSet,
)
from rest_framework.routers import DefaultRouter

from drf_spectacular.views import SpectacularAPIView, SpectacularSwaggerView

router = DefaultRouter()
router.register(r"samplecounts", SampleCountViewSet, basename="samplecounts")
router.register(r"metadata", MetadataViewSet, basename="metadata")
router.register(r"plates", PlateViewSet, basename="plates")
router.register(r"substrains", SubstrainViewSet, basename="substrains")
router.register(r"samples", SampleViewSet, basename="samples")
router.register(r"files", FileViewSet, basename="files")

urlpatterns = [
    path("admin/", admin.site.urls),
    path("api/auth/cookie/", CsrfCookieView.as_view(), name="auth-cookie"),
    path("api/auth/login/", LoginView.as_view(), name="login"),
    path("api/auth/logout/", LogoutView.as_view(), name="logout"),
    path("api/auth/users/me/", UserViewSet.as_view({"get": "retrieve"}), name="me"),
    path("api/download/<path:filepath>/", download_file, name="download_file"),
    path("api/", include(router.urls)),
    path("api/schema/", SpectacularAPIView.as_view(), name="api-schema"),
    path(
        "api/docs/",
        SpectacularSwaggerView.as_view(url_name="api-schema"),
        name="api-docs",
    ),
]
