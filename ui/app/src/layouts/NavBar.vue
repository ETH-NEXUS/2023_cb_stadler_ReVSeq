<script setup lang="ts">
import { ref } from 'vue';

const sideBarOpen = ref<boolean>(false);

const toggleSideBar = () => {
  sideBarOpen.value = !sideBarOpen.value;
};
</script>

<template>
  <section class="tw-relative tw-pb-24">
    <div
      class="tw-absolute tw-top-0 tw-right-0 tw-flex tw-w-full tw-h-3/4 md:tw-h-2/3 tw-bg-orange-200"
    >
      <img
        class="tw-absolute tw-bottom-0 tw-left-0 tw-right-0 tw-w-full tw--mb-4"
        src="../assets/background-orange-top.svg"
        alt="orange background"
      />
    </div>
    <nav
      class="tw-flex tw-mb-20 tw-justify-between tw-items-center tw-py-6 tw-px-10 tw-relative"
    >
      <a class="tw-text-lg tw-font-bold" href="#">
        <img
          class="tw-h-8"
          src="../assets/quasar-logo-vertical.svg"
          alt="logo image"
          width="auto"
        />
      </a>
      <div class="xl:tw-hidden">
        <button
          @click="toggleSideBar"
          class="tw-navbar-burger focus:tw-outline-none tw-text-indigo-900 hover:tw-text-indigo-800"
        >
          <svg
            class="tw-block tw-h-6 tw-w-6"
            fill="currentColor"
            viewbox="0 0 20 20"
            xmlns="http://www.w3.org/2000/svg"
          >
            <title>Mobile menu</title>
            <path d="M0 3h20v2H0V3zm0 6h20v2H0V9zm0 6h20v2H0v-2z"></path>
          </svg>
        </button>
      </div>
      <ul
        class="tw-hidden xl:tw-flex tw-absolute tw-top-1/2 tw-left-1/2 tw-transform tw--translate-x-1/2 tw--translate-y-1/2"
      >
        <li>
          <a
            class="tw-text-lg tw-mr-10 2xl:tw-mr-16 tw-font-extrabold hover:tw-text-indigo-800"
            href="#"
            >Product</a
          >
        </li>
        <li>
          <a
            class="tw-text-lg tw-mr-10 2xl:tw-mr-16 tw-font-extrabold hover:tw-text-indigo-800"
            href="#"
            >Features</a
          >
        </li>
        <li>
          <a
            class="tw-text-lg tw-mr-10 2xl:tw-mr-16 tw-font-extrabold hover:tw-text-indigo-800"
            href="#"
            >Pricing</a
          >
        </li>
        <li>
          <a
            class="tw-text-lg tw-font-extrabold hover:tw-text-indigo-800"
            href="#"
            >Resources</a
          >
        </li>
      </ul>
      <div class="tw-hidden xl:tw-flex tw-items-center">
        <a
          class="tw-inline-block tw-py-4 tw-px-6 tw-text-center tw-leading-6 tw-text-lg tw-text-white tw-font-extrabold tw-bg-indigo-800 hover:tw-bg-indigo-900 tw-border-3 tw-border-indigo-900 tw-shadow tw-rounded tw-transition tw-duration-200"
          href="#"
          >Log in</a
        >
      </div>
    </nav>
    <div class="tw-container tw-px-4 tw-mx-auto tw-relative">
      <div class="tw-max-w-5xl tw-mx-auto tw-text-center">
        <span
          class="tw-text-xl md:tw-text-2xl tw-font-extrabold tw-text-orange-500"
          >ETH Zürich</span
        >
        <h2
          class="tw-max-w-4xl tw-mx-auto tw-text-3xl sm:tw-text-4xl lg:tw-text-5xl tw-font-extrabold tw-font-heading tw-mt-1 tw-mb-6"
        >
          RevSeq
        </h2>
        <p
          class="tw-text-xl md:tw-text-2xl tw-font-extrabold tw-leading-8 tw-mb-12"
        >
          Gallia est omnis divisa in partes tres, quarum unam incolunt Belgae,
          aliam Aquitani, tertiam, qui ipsorum lingua Celtae, nostra Galli
          appellantur.
        </p>
        <div class="tw-flex tw-flex-wrap tw-mb-20 tw-justify-center">
          <a
            class="tw-inline-block tw-w-full md:tw-w-auto tw-mb-2 md:tw-mb-0 md:tw-mr-4 tw-py-4 tw-px-6 tw-text-center tw-leading-6 tw-text-lg tw-text-white tw-font-extrabold tw-bg-indigo-800 hover:tw-bg-indigo-900 tw-border-3 tw-border-indigo-900 tw-shadow tw-rounded tw-transition tw-duration-200"
            href="#"
            >Request a demo</a
          ><a
            class="tw-inline-block tw-w-full md:tw-w-auto tw-py-4 tw-px-6 tw-text-center tw-leading-6 tw-text-lg hover:tw-text-white tw-font-extrabold tw-bg-white tw-border-3 tw-border-indigo-900 hover:tw-bg-indigo-800 tw-shadow tw-rounded tw-transition tw-duration-200"
            href="#"
            >Learn More</a
          >
        </div>
      </div>
    </div>

    <div
      class="tw-navbar-menu tw-relative tw-z-50"
      :class="sideBarOpen ? '' : 'tw-hidden'"
    >
      <div
        class="tw-navbar-backdrop tw-fixed tw-inset-0 tw-bg-gray-800 tw-opacity-25"
      ></div>
      <nav
        class="tw-fixed tw-top-0 tw-left-0 tw-bottom-0 tw-flex tw-flex-col tw-w-full md:tw-w-5/6 tw-max-w-sm tw-py-8 tw-px-8 tw-bg-white tw-border-r tw-overflow-y-auto"
      >
        <div class="tw-flex tw-items-center tw-mb-8">
          <a
            class="tw-mr-auto tw-text-2xl tw-font-bold tw-leading-none"
            href="#"
          >
            <img
              class="tw-h-6"
              src="../assets/quasar-logo-vertical.svg"
              alt=""
              width="auto"
            />
          </a>
          <button class="tw-navbar-close" @click="toggleSideBar">
            <svg
              class="tw-h-6 tw-w-6 tw-text-gray-500 tw-cursor-pointer hover:tw-text-gray-500"
              xmlns="http://www.w3.org/2000/svg"
              fill="none"
              viewbox="0 0 24 24"
              stroke="currentColor"
            >
              <path
                stroke-linecap="round"
                stroke-linejoin="round"
                stroke-width="2"
                d="M6 18L18 6M6 6l12 12"
              ></path>
            </svg>
          </button>
        </div>
        <div class="tw-my-auto">
          <ul class="tw-py-10">
            <li class="tw-mb-1">
              <a
                class="tw-block tw-p-4 tw-text-lg tw-font-extrabold hover:tw-bg-gray-50 rounded"
                href="#"
                >Product</a
              >
            </li>
            <li class="tw-mb-1">
              <a
                class="tw-block tw-p-4 tw-text-lg tw-font-extrabold hover:tw-bg-gray-50 rounded"
                href="#"
                >Features</a
              >
            </li>
            <li class="tw-mb-1">
              <a
                class="tw-block tw-p-4 tw-text-lg tw-font-extrabold hover:tw-bg-gray-50 rounded"
                href="#"
                >Pricing</a
              >
            </li>
            <li class="tw-mb-1">
              <a
                class="tw-block tw-p-4 tw-text-lg tw-font-extrabold hover:tw-bg-gray-50 rounded"
                href="#"
                >Resources</a
              >
            </li>
          </ul>
        </div>
        <div>
          <a
            class="tw-block tw-py-4 tw-px-6 ttw-ext-center tw-leading-6 tw-text-lg tw-text-white tw-font-extrabold tw-bg-indigo-800 hover:tw-bg-indigo-900 tw-border-3 tw-border-indigo-900 tw-shadow tw-rounded tw-transition tw-duration-200"
            href="#"
            >Log In</a
          >
        </div>
      </nav>
    </div>
  </section>
</template>

<style scoped>
@import 'tailwindcss/base.css';
@import 'tailwindcss/components.css';
@import 'tailwindcss/utilities.css';
</style>
